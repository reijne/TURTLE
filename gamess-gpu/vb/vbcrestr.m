c
c  $Author: jvl $
c  $Date: 2012-03-25 14:32:50 +0200 (Sun, 25 Mar 2012) $
c  $Locker:  $
c  $Revision: 6255 $
c  $Source: /c/qcg/cvs/psh/GAMESS-UK/vb/vbcrestr.m,v $
c  $State: Exp $
c  $Log: vbcrestr.m,v $
c  Revision 1.40  2007/10/09 11:20:05  jvl
c  Minor updates follow:
c  - write(6,) => write(iwr,) to avoid printing outputs X (number of cores) times
c  - minor change to SUPER HYBRID (it has to follow HYBRIDS definition)
c  - truncated annoying warning messages to minimum
c  /marcin
c
c  Revision 1.39  2007/07/06 18:45:44  jvl
c  Minor bug fix in CRESTR with STRUC NORM.
c  Minor update of VBMO directive.
c  Addition of new stripblanks subroutine to strip all the front and trailing
c  blank spaces from string, while keeping all the within.
c  /marcin
c
c  Revision 1.38  2006/04/18 16:35:34  jvl
c  - added basis option to allow Vb to be run in anothr basis.
c  - cleaned up various bits of VB, like core partitioning
c  - adapted get1e
c
c  Revision 1.37  2006/01/13 17:56:48  jmht
c  Merging in changes from the release branch
c
c  Revision 1.36  2005/09/23 14:34:56  jmht
c  This checkin merges in the changes on the branch "release-7-0" back into the
c  development branch.
c
c  There are a number of changes to the way that the code is
c  configured. The configure script has been largely re-written and uses
c  the GNU config.guess script to determine the "GNU-triple" specifying
c  architecture of the machine. With names based on this there are then a number of
c  Machine-specific (MK) files stored in the config directory that hold
c  the compiler- and build-specific information for the different builds.
c  The MK files are inserted into Makefile.in in m4 to create the
c  final Makefile for the build. Makefile.in in m4 has also been
c  substantially changed to reflect this.
c
c  An couple of additional directories for validating the parallel code
c  have also been added. These are:
c
c  GAMESS-UK/examples/parallel_GAs and GAMESS-UK/examples/parallel_MPI
c
c
c  The large number of files in this checkin is largely due to the calls
c  to the intrinsic functions max0 and min0 being replaced with calls to
c  max and min.
c
c  The remaining larger changes are summarised below, with the log message for
c  the change being listed, together with the file where the change took
c  place (although obviously a single change can affect a number of files).
c
c  m4/analg.m:
c  Numerous changes to rationalise the current level of analysis and
c  printing when running "property atoms". This now includes generation of
c  properties in both global and principal axis, particularly for the
c  second moments and quadrupoles, plus modified field gradient O/P.
c  Results have been checked against the corresponding O/P from prop1e
c  (the code that is limited to s,p,d functions). - mfg
c
c  m4/basis.m:
c  We need to check on overflow of various arrays due to a job exceeding
c  the parameters in sizes BEFORE we call preharm. Otherwise the chances are
c  that the job dies with an obscure error message in preharm instead of
c  informing the user that her job exceeds the dimensions... - jvl
c
c  m4/c.m:
c  a significant re-working of this file in the hope that it will prove
c  a better starting point for moving around all the machines.
c  The use of memalign has been removed (as has the copy of the routine).
c  getmem will check if alignment is good and will exit(-1) with a
c  message if not, so we should soon find if it memalign was actually
c  needed
c  Declarations of malloc have been removed and a more standard set
c  of #includes introduced. It should be easier to maintain these.
c  Some ipsc and dec specific code has been removed -psh
c
c  m4/cphf.m:
c  Remco had discovered a bug that caused calculations with
c     dft
c     runtype hessian
c  and
c     runtype hessian
c     dft
c  to give different results. This has been fixed now by:
c  a) controling the selection of the CPHF solver in a different way, before
c     the logical scalar was used, now CD_active() is used directly
c  b) the invokation of CD_gradquad(.true.) for hessian calculations has been
c     put in a different location to ensure it is always invoked for hessian
c     jobs.
c
c  - hvd
c
c
c  m4/ga.m:
c  Fixed bug 23: The timing errors were caused by the subroutines wrt3_ga
c  and rdedx_ga exiting without call end_time_period. This happened in
c  the "special case for ed19 c/z vectors". Fixed this by replacing the
c  return statement with goto 999, and adding a corresponding continue
c  statement just before the end_time_period calls. - hvd
c
c  This is a fix for the bug that was causing problems when restarting ga jobs from
c  a dumpfile on disk. For zero word writes, a return statement had been added in wrt3_ga
c  and rdedx_ga as this was needed for ed16 under the CI code. Unfortunately this meant that the
c  block counter wasn't updated for zero length blocks when reading in a dumpfile, causing all
c  blocks following an empty block to be displaced by one. The fix is just to make the return
c  statement for zero word writes conditional to ed19. - jmht
c
c
c  m4/guess.m:
c  Got the *definite* version of atoms; the charge is concentrated in the
c  highest (open) shell, as described in the paper -jvl
c
c  corrected nattempt to uswe charge and spin simulaneously in atdens.
c  is not needed and now forbidden - jvl
c
c  The message "error in atom scf" is most often triggered by forgetting
c  to specify an ECP with an ECP basis set, or by making a mistake in the
c  basis set input. This is not very clear from the error message
c  however, therefore I have added some text to guide the user in finding
c  the likely cause of the error. - hvd
c
c  m4/index.m:
c  replace the explicit zora reads by zora calls in index4 and tran4
c  and restricted usage to (reommended) atomic zora
c  also added an sgmata_255 in index4 allthough th resulting fock matrix
c  seems to be ignored - jvl
c
c  m4/integb_lib.m:
c  An attempt to trap faulty ECP directives. This traps sillies like
c      PSEUDO STRSC
c  which does not mean anything to the code, but upto now the STRSC string
c  was silently ignored. This would confuse users as they would think they had
c  clearly specified which ECP they wanted. Now the code prints the acceptable
c  forms of the initial ECP card and stops. Hopefully this will trigger the
c  inexperienced user to read the manual (one can always hope). -hvd
c
c  machscf.m:
c  disable global array file system by default on all platforms.
c  there are still bugs in this code so we are disabling it until
c  it is fixed and then we'll discuss if it should be re-enabled
c  for specific platforms. - psh
c
c  mains.m:
c  added a list of M4 keywords to the build information that is
c  incorporated in the executable.
c  This is not printed by default (unlike the date/user of the build)
c  but can be obtained by including the "keys" directive. The idea is
c  that it should help us work out exactly which options are active
c  for a particular build. - psh
c
c  A number of changes were made to the calculations of frequencies
c  to achieve more consistent results among different parts of the code.
c  The main changes are:
c  1) Initialise the physical constants in common/funcon/ from the values
c     stored in common/phycon/ to ensure better consistency.
c     (subroutine pacini [master.m])
c  2) Calculate conversion factors from physical constants where needed.
c     (subroutine iranal [anale.m], subroutine rotcon [optim.m],
c      subroutine dr2fre [sec2e.m])
c  3) Use repeated modified Gramm-Schmidt orthonormalisation instead of
c     less accurate modified Gramm-Schmidt or even Gramm-Schmidt.
c     (subroutine prjfco [anale.m], subroutine vibfrq [optim.m])
c  4) Applying the projection to remove translational and rotational
c     coordinates to the mass weighted force constant matrix instead
c     of the normal force constant matrix.
c     (subroutine prjfcm [anale.m])
c  5) Modified the subroutine mofi [optim.m] so that it now returns
c     the centre-of-mass, the moments of inertia, and the eigenvectors
c     of the moments of inertia tensor.
c  6) Use subroutine mofi to compute the moments of inertia where needed.
c     (subroutine vibfrq [optim.m], subroutine thermo [optim.m]) - hvd
c
c
c  m4/newutil3:
c  newutil3 option (ie selecting additional blas-like routines by individual
c  names rather than based on the machine) is now the adopted mechanism.
c  The code corresponding to the old way of doing it has been removed. - psh
c
c  Revision 1.35.2.1  2005/07/19 06:52:19  wab
c
c  min0/max0 changes (see m4 checkin)
c
c  Revision 1.35  2005/05/11 22:37:31  jvl
c  made paralleol vb - I8 completer + few checks and better errormessages
c  sendrecv introduced to avoid race conditions (added to ga-stuff really)
c
c  Revision 1.34  2005/04/22 11:07:52  jvl
c  All vb criteria are now in common vbcri and may be set
c  in the input. Also their defaults are now consistent
c
c  Revision 1.33  2005/03/25 23:37:08  jvl
c  removed i8 error message
c
c  Revision 1.32  2005/02/05 18:04:25  wab
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
c  Revision 1.31  2004/03/27 15:04:35  jvl
c  corrected multi-structure struc + wrong zero call
c
c  Revision 1.30  2004/03/24 17:26:07  jvl
c  ome fortran corrections .....
c
c  Revision 1.29  2004/03/24 14:57:45  jvl
c  fixed a bug in print-out of branching diagrams and leading terms (JJE).
c
c  Revision 1.28  2004/03/23 15:07:53  jvl
c  enhancved conf manual and added default sections
c
c  Revision 1.27  2004/03/15 12:31:39  jvl
c  Made conf manual (or better) styruc so to accept more structures
c  Now you really can do a garbge calculation......
c
c  Revision 1.26  2004/03/13 22:57:49  jvl
c  reinstated the noproj option in crestr
c
c  Revision 1.25  2003/10/08 17:31:10  jvl
c  several bugfixes (JJE 8/10/2003).
c
c  Revision 1.24  2003/08/27 10:40:33  jvl
c  bugfix in prweig (JJE 27/8/2003)
c
c  Revision 1.23  2003/08/26 13:35:56  jvl
c  Weights are printed together with their branching or rumer diagram.
c  (JJE 26/8/2003)
c
c  Revision 1.22  2003/08/22 14:40:58  jvl
c  Fixed bug in leading term print. Leading terms and bonding arrangements
c  are stored now on ed7 (JJE 22/8/2003).
c
c  Revision 1.21  2003/04/07 15:35:11  hvd
c  Zero-length strings are a fatal error with the WATCOM compiler. Fixed this
c  in subroutine spinef.
c
c  Revision 1.20  2003/04/07 14:47:10  hvd
c  'conf man' led to a goto into the 'struc' branch of the directive
c  if-then-elseif-statement. this is a fatal error when using the
c  WATCOM compilers. The if-statement has been revised to remove the
c  goto-statement.
c
c  Revision 1.19  2003/03/23 16:30:21  jvl
c  fixed omitted initialisation (relevant for determinant printing)
c  official error message for spin off  ; ytest => itest in servec
c
c  Revision 1.18  2003/02/18 17:17:13  jvl
c  A lot of changes, in general:
c  - In vbci the call to bsmat for natural orbitals is no longer necessary
c  - The suboptions pert and fock for optimise in scf are implemented
c  - The subroutine vbscf is completely restructured and is more logical
c  - Addition of the makepsi0, which uses a small extra 2-electron,
c    transformation
c    JJE - 18 feb 2003
c
c  Revision 1.17  2002/09/05 14:24:28  jvl
c  Print of leading term modified. When nsel>40 no printing of leading term,
c  only a warning (JJE).
c
c  Revision 1.16  2002/05/28 15:07:47  jvl
c  General cleanup
c  got rid of all fortran-files (now part of ed7 is used)
c  handle 3 more common's through include
c  make vector sections in vb more logical
c  added print parameter to getqvb
c  ....
c
c  Revision 1.15  2002/03/18 17:09:38  jvl
c  Ridiculous bugfix (JJE)
c
c  Revision 1.14  2002/03/04 12:07:39  jvl
c  In CRESTR when printflag is "m" or "all" the leading term is printed
c  behind the structure number and the number of determinants (JJE).
c
c  Revision 1.13  2001/11/22 15:56:19  jvl
c  corrected determinant wiping (did not work for branching diagrams)
c  and changed isign usage and arrays which were dynamic
c
c  Revision 1.12  2001/10/24 21:45:38  jvl
c  corrected few problems in geometry optimisation with VB and added initial stuff
c  to limit output later..
c
c  Revision 1.11  2001/10/19 16:30:22  jvl
c  clean hybrid (default) now tries to make hybrid exact
c  root broadcasts orbitals + excitation info for parallel jobs
c  hybrids are checked for consistency
c  infato now contains fragment/ao (atomao) numbers making lots of routines easier
c  some routines had name-changes (clao=> clvecvb)
c  eliminated confusing excitation prints
c
c  Revision 1.10  2001/10/15 15:32:39  jvl
c  major parallel overhaul; makes asynchronous switchable
c  no luck yet however ...........
c
c  Revision 1.9  2001/09/19 16:02:42  jvl
c  some corrections resulting from DEC compile
c
c  Revision 1.8  2001/08/25 22:46:20  jvl
c  max => maxi_crestr (clashed with peigss)
c  remove dependency on mpistaus etc from vbtran
c
c  Revision 1.7  2001/07/04 21:38:07  jvl
c  few more irritatingf prints corrected
c
c  Revision 1.6  2001/06/27 15:57:36  jvl
c  bugfix in subroutine remdts. Before this routine only worked when 1 confi-
c  guration was specified. Fixed (JJE).
c
c  Revision 1.5  2001/06/12 12:21:16  jvl
c  - Removed several bugs
c  - Improved parallel 4-index transformation including changes in getin2
c  - Added timing for VB routines
c
c  Revision 1.4  2001/04/14 19:40:07  jvl
c
c  changes b from linuxppc port + some vb fixes
c
c  Revision 1.3  2001/02/08 14:34:35  jvl
c  Adapted parallelisation for GAMESS-UK.
c  An MPI and GA version have been created
c
c  Revision 1.2  2000/03/31 14:18:25  jvl
c  Changed vnorm to vnrm
c  Changed get1e to use getmat
c  Changed scopy to icopy
c  Improved on information about titles
c  Fokke Dijkstra
c
c  Revision 1.1  2000/02/16 12:20:42  jvl
c  adding vb files to gamess repository
c
c Revision 1.5  1997/05/30  11:10:21  fokke
c Changed icopy to icopy and removed variables icopy
c
c  Revision 1.9  1999/06/03 12:33:12  fokke
c  subroutine redund has been removed
c
c  Revision 1.8  1999/05/21 16:37:25  joop
c  latest fixes by arno (klogic / print)
c
c  Revision 1.7  1999/05/19 13:47:20  joop
c  msocc (# of singly occupied orbitals) parameter enhanced
c  Arno Blok
c
c  Revision 1.6  1999/05/17 13:33:45  joop
c  changed determninant generating to only generate determinant as they are needed.
c  the logic numbers (original order number) is now stored in seperate array klogic
c  In this array the determinants are now found by logic number
c  This saves lots of time and memory for bigger cases.
c  *note if all rumer (or branching) diagram (or the last) are required still
c  everything is generated. For Rumer this is fixable ...........
c  Arno Blok ... May 1999
c
c  Revision 1.5  1997/05/30 11:10:21  fokke
c  Changed icopy to scopy and removed variables icopy
c
c Revision 1.4  1997/05/22  12:53:09  joop
c changed imove to scopy
c
c Revision 1.4  1997/05/22  12:53:09  joop
c changed imove to scopy
c
c Revision 1.3  1997/05/22  11:31:48  joop
c Added include macro.cpp
c
c Revision 1.2  1996/10/29  15:38:29  joop
c rcs info + structures file + dimensions extended
c
c
      subroutine addcol (karr1,ndima1,ncola1,karr2,ndima2,ncola2,nrowa2)
c***********************************************************************
c
      implicit REAL  (a-h,o-z) , integer (i-n)
c
c***********************************************************************
c
c     'addcol' adds all the columns of array 'karr2' to array 'karr1'.
c
      dimension karr1 (ndima1, * ),
     &          karr2 (ndima2, * )
c
      do 10 iexcol=1,ncola2
         do 20 iexrow=1,nrowa2
            karr1(iexrow,ncola1+iexcol) = karr2(iexrow,iexcol)
20       continue
10    continue
c
      return
      end
      subroutine chkcon(kconf,ndim,nelec,nconf,ndelet,kdelet,nspath)
c***********************************************************************
c
      implicit REAL  (a-h,o-z) , integer (i-n)
c
c***********************************************************************
c
c     'delcol' may delete some columns of the array 'kconf', so that
c     only unique configurations remain. row-element 'nelec+1' of each
c     column equals the number of double occupied orbitals in the cor-
c     responding configuration.
c     the number of deleted columns is given through to the main program
c     by means of dummy argument 'ndelet'.
c
      dimension kconf  ( ndim , * ),
     &          kdelet ( nconf ),
     &          nspath ( * )
c
      jcount = 0
      do 1 i=1,nconf
         kdelet(i) = 0
1     continue
c
      do 10 icol1=1,nconf-1
         do 20 icol2=icol1+1,nconf
c
            do 30 irow=1,nelec
               if (kconf(irow,icol1).ne.kconf(irow,icol2)) goto 40
30          continue
            if (nspath(icol1).ne.0.or.nspath(icol2).ne.0) goto 40
c
c     column 'icol1' is equal to column 'icol2', check wether column
c     'icol2' hasn't been registrated before
c
             do 50 idel=1,jcount
                if (icol2.eq.kdelet(idel)) goto 40
50           continue
c
c     registrate the twice occuring column (only one of them must be
c                                                            deleted)
c
                    jcount = jcount + 1
            kdelet(jcount) = icol2
c
40          continue
20       continue
10    continue
      ndelet = jcount
c
      return
      end
      subroutine colint(kconf,melec1,nrefco,kinner,ninner)
c...
c...  store internal orbitals, make sure kinner(i+1) > kinner(i) holds
c...
      implicit REAL (a-h,o-z), integer (i-n)
c
      common /indep/ time0,nelec,mult
c
      dimension kconf (melec1,*),
     &          kinner(*)
c
c...  set up ordered list of internals from the first configuration
      ndocc = kconf(nelec+1,1)
      nsocc = nelec - 2*ndocc
      ninner = 0
      do 10 iorb=1,nelec
c...  (store double occupied orbitals only once)
         if (iorb.gt.nsocc.and.((iorb+1)/2)*2.ne.iorb) go to 10
         kinner(ninner+1) = kconf(iorb,1)
         ninner = ninner + 1
10    continue
      call ordernumbs(kinner,ninner)
c...  add unregistrated orbitals of remaining configurations to the list
      do 20 iconf=2,nrefco
         ndocc = kconf(nelec+1,iconf)
         nsocc = nelec - 2*ndocc
         do 30 iorb=1,nelec
            if (iorb.gt.nsocc.and.((iorb+1)/2)*2.ne.iorb) go to 30
c...  find a position for the orbital on the list, that maintains the 
c...                                                            ordering
            jpos = 1
40          if (jpos.le.ninner) then
               if (kconf(iorb,iconf).eq.kinner(jpos)) then
c...  the orbital is on the list already
                  go to 30
               else if (kconf(iorb,iconf).gt.kinner(jpos)) then
                  jpos = jpos + 1
                  go to 40
               end if 
            end if
c...  add orbital to the list
            call open1(kinner,jpos,ninner+1)
            kinner(jpos) = kconf(iorb,iconf)
            ninner = ninner + 1
30       continue
20    continue
c
      return
      end
      subroutine copi(ifrom,ito,nelem,incr)
c...
c...  copy integer arrays. the copy-order is reversed by a negative 
c...  increment. (proper choice of ifrom, ito and incr ensures
c...  proper copying if ifrom and ito belong to the same array)
c...  ifrom and ito must be called with their start adresses
c...
      implicit REAL (a-h,o-z), integer (i-n)
c
INCLUDE(../m4/common/iofile)
c
      dimension ifrom(*), ito(*)
c
      if (incr.lt.0) then
         jfirst = nelem
         jlast  = 1
      else
         jfirst = 1
         jlast  = nelem
      end if
      do 10 icp=jfirst,jlast,incr
         ito(icp) = ifrom(icp)
10    continue
c
      return
      end

      subroutine copl(lfrom,lto,nelem,incr)
c...
c...  copy logical arrays. the copy-order is reversed by a negative 
c...  increment. (proper choice of lfrom, lto and incr ensures
c...  proper copying if lfrom and lto belong to the same array)
c...  lfrom and lto must be called with their start adresses
c...
      implicit REAL (a-h,o-z), integer (i-n)
c
      logical lfrom(*), lto(*)
c
      if (incr.lt.0) then
         jfirst = nelem
         jlast  = 1
      else
         jfirst = 1
         jlast  = nelem
      end if
      do 10 icp=jfirst,jlast,incr
         lto(icp) = lfrom(icp)
10    continue
c
      return
      end

      subroutine copr(rfrom,rto,nelem,incr)
c...
c...  copy real arrays. the copy-order is reversed by a negative 
c...  increment. (proper choice of rfrom, rto and incr ensures
c...  proper copying if rfrom and rto belong to the same array)
c...  rfrom and rto must be called with their start adresses
c...
      implicit REAL (a-h,o-z), integer (i-n)
c
      dimension rfrom(*), rto(*)
c
      if (incr.lt.0) then
         jfirst = nelem
         jlast  = 1
      else
         jfirst = 1
         jlast  = nelem
      end if
      do 10 icp=jfirst,jlast,incr
         rto(icp) = rfrom(icp)
10    continue
c
      return
      end
      subroutine crestd(kconf,knstru,kntdet,kstruc,
     *                  coeff,kdeter,klogic,nspath,ispath,maxsp,
     *                  coefff,kdetps,idetps,pacdet,
     *                  kscr  ,mscr  ,scr   ,nelec1,nconf ,ntdet,
     *                  irumer)
c
      implicit REAL (a-h,o-z), integer (i-n)
c
c...    driver and input-processor for crestr
c
INCLUDE(common/turtleparam)
INCLUDE(common/c8_16vb)
c
      logical    symc
c
      character *4  ied
      character *8  acct, date, time, jobnam
      character *8  label
      character *15 object
      character *8  group, irrep
c
      common/crestn/ struct
      character*44 struct
c
      dimension pacdet (*),
     &          kconf  (nelec1,nconf),
     &          knstru (nconf),
     &          kntdet (nconf),
     &          coefff (*),
     &          kdetps (*),
     &          idetps (*)
      dimension kstruc (*),
     &          coeff  (*),
     &          kdeter (nelec1-1,*),
     &          klogic (*),
     &          nspath (nconf),
     &          ispath (maxsp,nconf),
     &          kscr   (mscr),
     &          scr    (*),
     &          irumer (*),
     &          idoub(mxorbvb),
     &          iarr(mxorbvb)
c
      common /indep/  time0,nelec,mult
INCLUDE(common/logics)
      common /hans/   msocc,mstrco,mdetco,mdetst,mperm,mstruc
INCLUDE(../m4/common/iofile)
INCLUDE(../m4/common/timeperiods)
INCLUDE(common/aivb)
c       
      if(manual) then
c...
c...  there  are just 'manually defined' structure
c...  note right order is  aaaaabbbb
c...
c...  set ndoub and idoub (doub. indexes) values
        do i=1,mxorbvb
           iarr(i)=0
        enddo
c...  check how many times an orbital is mentioned in the complete set
        iind = 1
        do i=1,nconf
          do j=1,nelec
            iarr(kdeter(j,iind))=iarr(kdeter(j,iind))+1
          enddo
          iind = iind + nspath(i)
        enddo
c...  if an orbital is mentioned 2 times the number of configurations
c...  it is doubly occupied in all configurations (so nconf*2).
        ndoub=0
        do i=1,mxorbvb
           if (iarr(i).eq.nconf*2) then
              ndoub=ndoub+1
              idoub(ndoub)=i
           endif
        enddo

         nalpha = (nelec + mult-1)/2
         ntstru = nconf
         ncfff  = ntdet
         jpacd = 1
         jkdeter = 1
         is = 1
         maxdet = 0
         maxstr = 0
         do i=1, nconf
            knstru(i) = 1
            kntdet(i) = nspath(i)
            kdetps(i) = nspath(i)        
            call pack(pacdet(jpacd),n8_16,kdeter(1,jkdeter),
     1                nelec*kntdet(i))     
            jpacd = jpacd + (nelec*kntdet(i)-1)/(64/n8_16) + 1
            jkdeter =jkdeter + kntdet(i)
            do j=1,kntdet(i)
               idetps(is) = j
               is = is + 1
            end do
         end do
         npd = jpacd - 1
c                                               
         imax   = 0     
         do j=1,ntdet
            do i=1,nelec
               imax = max(kdeter(i,j),imax)
            end do                               
         end do
         maxdet = ntdet
         maxstr = nconf
c     
c...     output
c
         call strtrm(struct,nn)
         write(iwr,215) nconf,ntstru,ntdet,struct(1:nn)
c
         jpacd = 1
         if ( aivb_set .eqv. .true. ) then
           write(iwr,'(a)') 'Automatically generated structures:'
         else
           write(iwr,'(a)') ' Manually supplied structures:'
         end if
         do i=1,nconf
            write(iwr,'(a,i3,32a)') ' == structure',i,' ==',
     1      (' a ',l=1,nalpha),(' b ',l=nalpha+1,nelec)
             if ( aivb_set .eqv. .true. ) then
               write(iwr,'(a,a,a)') ' Atom ',
     1         'configurations: ',aivb_confdescr(i,1)
               write(iwr,'(a,a)') ' Atomic states: ',
     1         aivb_confdescr(i,3)
             end if
            do j=jpacd,kntdet(i)+jpacd-1 
               write(iwr,'(a,f10.5,2x,(20i3))') ' coef ',coefff(j),
     1                                         (kdeter(k,j),k=1,nelec)
            end do
            if ( i .lt. nconf ) write(iwr,'(a)') ' '
            jpacd = jpacd + kntdet(i)
         end do
c
      else     
c
c                      ***** configuration-cycle *****
c
         stime0 = cpulft(1)
         ntdet  = 0
         ntstru = 0
         kkd    = 1
         jpacd  = 1
         ks     = 1
c
         do 40 iconf=1,nconf
c
c               create structures (based on leading terms)
c
            if (allout) then
               write(iwr,3) iconf
3              format (/,' ',11x,'**** build structures from ',
     &                 'configuration',i4,' ****')
            endif     
c
c...     set highest spin path
c
            msp = 0
            do isp=1,maxsp
               msp = max(ispath(isp,iconf),msp)
            end do
c   
            call spinef (kconf(1,iconf),nconf,kdeter,klogic,
     &                   kntdet(iconf),kstruc,
     &                   coeff,kndets,knstru(iconf),allout,projec,
     &                   kscr,mscr,mdetst,nelec,msp,irumer)
            if (knstru(iconf).eq.0)  
     &                      call vberr(' spinef rejected configuration')
c
c...     produce final results  for structures
c
                       call schmic(coeff,mdetst,kndets,knstru(iconf),
     *                             kstruc,scr,kntdet(iconf),
     *                             coefff(kkd),kdetps(ntstru+1),
     *                             idetps(kkd),ncc,
     *                             nspath(iconf),ispath(1,iconf),
     *                             rumer(iconf),kdeter,nelec,klogic,
     *                             irumer)

cmarcin      write(iwr,'(2(a,I3))')
cmarcin     1      ' nnwdet/kntdet: ',kntdet(iconf),
cmarcin     2      '  nstruc/knstru: ',knstru(iconf)
cmarcin      write(iwr,'(a,20I3)') ' kdetps: ',(kdetps(i), i=1,mstruc)
cmarcin      write(iwr,'(a,20I3)') 
cmarcin     1          ' idetps: ',(idetps(i), i=1,ncc)
cmarcin      write(iwr,'(a,90F8.4)')
cmarcin     1          ' coefff: ',(coefff(i), i=1,ncc)
c
c...     remove unused determinants
c
            if (kntdet(iconf)*2.gt.mscr) call vberr('remdt memory')
            call remdts(kdeter,nelec,kntdet(iconf),idetps(kkd),
     *                  kdetps(ks),knstru(iconf),
     *                  kscr,kscr(kntdet(iconf)+1))
            ks = ks + knstru(iconf)
c
c...     pack determinants
c
            call pack(pacdet(jpacd),n8_16,kdeter,nelec*kntdet(iconf))     
            jpacd = jpacd + (nelec*kntdet(iconf)-1)/(64/n8_16) + 1
            ntdet =  ntdet + kntdet(iconf)
            ntstru = ntstru + knstru(iconf)
            kkd = kkd + ncc
c
         if (ntstru.gt.maxcon) call dimerr(' vb structures',maxcon)

40       continue
         write(iwr,55) cpulft(1) - stime0
55       format(/,t10,'** generation of spin functions in',f8.2,
     &                                        ' cpu-seconds **') 
         ncfff = kkd - 1
         npd = jpacd - 1
c
c                 ****  end  configuration-cycle *****
c
c       final output is
c       pacdet(jpacd-1) : gepackte determinanten
c       kntdet(nconf)   :  determinants / config
c       knstru(nconf)   :  structures / config
c       kdetps(ntstru)  :  determinants / structure
c       idetps(kkd-1)   : indices of determinants for a structure
c       coefff(kkd-1)   : coeff's of determinants for a structure
c
c          some arrangements, necessary if spin-projection is off
c
        if (.not.projec) then
            is = 1
            do 50 i=1,nconf
               knstru(i) = kntdet(i)
               do 45 j=1,kntdet(i)
                  kdetps(is) = 1
                  idetps(is) = j
                  coefff(is) = 1.0d0
45             is = is + 1
50          continue
            ntstru = is - 1
            ncfff = is - 1
        endif
c
         imax = 0
         maxdet = 0
         maxstr = 0
         do 19 i=1,nconf
            
            maxdet = max(kntdet(i),maxdet)
            maxstr = max(knstru(i),maxstr)
            do 21 j=1,nelec
              imax = max(kconf(j,i),imax)
21          continue
19       continue
c
         call strtrm(struct,nn)
         write(iwr,215) nconf,ntstru,ntdet,struct(1:nn)
215      format(///1x,55(1h*),/,t11,'final results',/,
     *           t10,'number of configurations : ',i6,/
     *           t10,'number of structures     : ',i6,/
     *           t10,'number of determinants   : ',i6,/
     *           t10,'written to the file      : ',a,/
     *           1x,55(1h*))
c
         if (allout.or.limout) then
c
c               write all structures and determinants to tape6 (output)
c
            ks = 1
            kp = 1
            kc = 1
            write(iwr,65)
65          format(/' == list of configurations == ')
            do 340 iconf=1,nconf
c
               if((nelec*knstru(iconf)).gt.mscr) call corfait((nelec*
     &          knstru(iconf)),mscr,'not enough memory before prconf')
               call leadsup(irumer((ks-1)*nelec+1),kscr,
     &          (nelec*knstru(iconf)))
               call prconf(knstru(iconf),kntdet(iconf),pacdet(kp),
     *                     kdetps(ks),idetps(kc),coefff(kc),nc,
     *                     kdeter,nelec,kconf(1,iconf),iconf,
     *                     kscr)
               ks = ks + knstru(iconf)
               kc = kc + nc
               kp = kp + (nelec*kntdet(iconf)-1)/(64/n8_16) + 1
c
340         continue
c
            write(iwr,68)
68          format(/' ======== end of list  ========= ')
         end if
c
c
      end if
c
c...  store bonding patterns on ed7
c
      k7bond = kscra7vb('k7bond',ntstru*nelec,'i','w')
      call wrt3is(irumer,ntstru*nelec,num8)
c
c              write all relevant information to tape15
c
      call tape15(nelec,ntdet,ntstru,nconf,kconf,nelec+1,
     &            pacdet,kntdet,knstru,kdetps,idetps,coefff,ncfff,
     &            npd,mult,imax,maxdet,maxstr,ndoub,idoub)
c
      write(iwr,66) cpulft(1) - time0
      call end_time_period(TP_VB_STRUC)
66    format(/,t10,'** total cpu time used: ',f8.2,' cpu-seconds **')

c
      return
      end
      subroutine cresti(kconf ,melec1,nspath,ispath,maxsp ,korb  ,
     &                  kvirt ,idconf,krejco,kscr  ,mconf ,nconf ,
     &                  nsp0  ,ninner,nvirt ,coefff,kdeter,mdet  , 
     &                  ntdet ,v, qq )
c...
c...     input-processor for crestr     
c...
c
c                                input 
c
c     all directives start with a key-word, often followed by one or 
c     more corresponding input-lines. the order in which the directives
c     appear is arbitrary, but the obligatory nelec-directive must be  
c     first and the enter-directive must be last.
c  
c     orbitals are represented by integers. a double occupied orbital
c     is defined by the double occurance of an integer. a string of
c     numbers may or may not be closed by the closing tag 'end'. 
c     integer-strings may always be spread out over several lines, 
c     except where indicated.
c
c                              directives
c
c---- the print-directive
c     this directive consists of the key-word 'print', followed by an 
c     output-mode specifier. if this is 'all', then a large amount of
c     output, considering the generation of structures, is given. if it 
c     is 'limi' or 'm', just the defined structures and determinants are
c     shown. default only very modest information is given.
c
c---- the multiplicity-directive
c     input-lines
c                  - 'mult'
c                  - the multiplicity (2s+1)
c     if ommited, the multiplicity is set equal to 1 for equal number of
c     electrons and equal to 2 for odd number of electrons
c
c---- the symmetry-directive
c     input-lines
c                  - 'sym'
c                  - the number of basis functions
c                  - the pointgroup-symbol, the irrep.-symbol of the
c                                                             molecule
c                  - an irrep.-symbol, the corresponding basis functions
c                    (more of such lines can occur, until the irreps of
c                     all relevant basis functions have been established
c     if ommited, no use can be made of symmetry.
c
c---- the spin-directive
c     this directive consists of the key-word 'spin', followed by a 
c     spin-projection specifier. the options 'rumer' and 'lt' give rise
c     to rumer-functions, the options 'bd', 'yk' and 'gen' to branching
c     diagram-functions. if the option 'off' is set, the determinants
c     will be defined as independent entities, and will not be gathered
c     in configurations or structures.
c     default rumer-functions will be defined
c
c---- the configuration-directive
c     with this directive ready-made configurations can be given.
c     input-lines
c                  - 'conf'
c                  - a configuration
c                    (more of such lines can occur)
c                  - 'end'
c     the definition of one configuration may be spread out over several
c     lines, but a new configuration must always appear on a new line
c     any configuration may be followed by a spin-path specifier. 
c     default all structures corresponding to a configuraton are formed,
c     by performing the projection, set with the spin-directive. if a 
c     spin-path specifier follows a configuration, the normal projection 
c     is overruled for the particular configuration. the same options 
c     may be used as with the spin-directive (except for the option 
c     'off'). if the option is followed by one or more integers only the
c     the selected spin-paths will be defined. the selection-request
c     must be defined entirely on the same line as the last orbital of
c     the configuration
c     example:
c              spin bd
c              conf
c              1 1 2 3 (end)
c              1 1 2 4 (end) rumer (end)
c              1 1 2 5 (end) bd 2 (end)
c              end  
c     all bd-functions will be made for the first configuration, and all 
c     rumer-functions for the second configuration. for the third confi-  
c     guration only the second bd-function will be made.
c
c---- the configuration-directive full-ci extenstion
c     with this directive all possible combinations of 'm' electrons 
c     distributed over 'n' orbitals, will be generated
c     input-lines
c           - 'conf'
c           -  'all' 'number of orbitals' 'number of electrons'
c           - 'end'
c     sometimes, in case of molecules with more than one hydrogen,
c     one would like to restrict the configurations to those, where
c     no more that one hydrogen at once is doubly occupied
c     this can be done by adding
c           - 'restr' 'num. of orb. to be restrict' 'orb. num. to restrict'
c     in the same line as 'all', after 'number of electrons'
c
c---- the man-subdirective
c.... (might change in future)
c     with this directive structures can be given by hand
c     input-lines
c                  - 'stru' 'conf man'
c                  - coefficient determinant (more of such lines may
c                                                             occur)
c                  - end
c     example: (singlet)
c              stru
c              1.0 1 2 1 3
c              1.0 1 3 1 2
c              end
c     the combined used of the man-directive with the conf- is
c     strongly dissuated
c
c---- the mrsd-directive
c     with this directive a group of reference configurations and a set
c     of virtuals can be defined, with which single and double excita-
c     ions will be performed
c     input-lines:
c                  - 'mrsd' 'stru'/'conf' 
c                    'refer'
c                  - a reference configuration (followed by 'end')
c                    (more of such lines may occur)
c                  - 'virt'
c                  - a set of virtuals (followed by 'end')
c                  - 'end'
c     the mrsd-directive cannot be used in combination with other
c     configuration-generating directives. the order of refer- and virt-
c     directive is arbitrary
c
c---- the enter-directive
c     this obligatory directive must always end an input-section and 
c     consists of the single word 'enter'
c   
      implicit REAL (a-h,o-z), integer (i-n)
c
INCLUDE(common/splice)
INCLUDE(common/vbcri)
INCLUDE(common/c8_16vb)
INCLUDE(common/turtleparam)
c
      logical symc
      logical spin
      logical messag
      logical next,mfirst,onorm
c
      character*8  label, label1
      character*15 object
c92      character*8  group,irrep
      character*4  group,irrep
c
      dimension kconf (melec1,mconf),
     &          krejco(mconf),
     &          ispath(maxsp ,mconf),
     &          nspath(mconf)
      dimension kvirt (*),
     &          coefff(*),         
     &          kdeter(melec1-1,*),
     &          korb(*),
     &          idconf(4,*),
     &          kscr  (*),
     &          qq(*),
     &          v(*)
c
      common /indep/  time0,nelec,mult
      common /hans/   msocc,mstrco,mdetco,mdetst,mperm,mstruc
INCLUDE(common/logics)
INCLUDE(../m4/common/iofile)
INCLUDE(common/allconf)
c
INCLUDE(common/aivb)
      external aivbdata
c
      character*1 figure(10)  
c
      data (figure(i), i=1,10)   
     &     /'1','2','3','4','5','6','7','8','9','0'/     
c
33    format(//30x,a)
c...
c...  initialize 
c...
      allout = .false.
      limout = .false.
cjje  mult   = 0
      symc   = .false.
      projec = .true.
      spin   = .true.
      nconf  = 0
      manual = .false.
      mrsd   = .false.
      mfirst = .true.
      cridep = 1.0d-13
      do i=1,mconf
         do j=1,maxsp
            ispath(j,i) = 0
         end do
      end do
c
c...
c...  *********************** input-cycle *****************************
c...
1     if (mfirst) then
         mfirst=.false.
      else
         call input
      end if    
      call inpa(label)
2     continue
c
      if (label(1:4).eq.'255+'.or.label(1:4).eq.'>255') then
c...
c...  allow more then 255 orbitals - use 16 bit packing
c
         n8_16 = 16
c
      else if (label(1:4).eq.'255-'.or.label(1:4).eq.'<255') then
c...
c...  allow less then 255 orbitals - use 8 bit packing
c...
         n8_16 = 8
c
c
      else if (label.eq.'print') then
c...
c...                        set output-mode
c...
         call inpa(label)
         if (label(1:3).eq. 'all'.or.label(1:1).eq.'h') then
            allout = .true.   
            limout = .true.
         else if (label(1:4).eq.'limi'.or.label(1:1).eq.'m') then
            limout = .true.
         else
            write(iwr,33) label
            call vberr(' unrecognized output-mode')   
         end if
c
C
C*** TO BE SUPPLIED BY GAMESS
C
c        else if (label(1:4).eq.'mult') then
c...
c...                       read multiplicity
c...
c        call inpi(mult)
c        if (mult.le.0) call vberr(' impossible multiplicity')
c        if (((mult/2)*2.eq.mult.and.(nelec/2)*2.eq.nelec).or.
c    &      ((mult/2)*2.ne.mult.and.(nelec/2)*2.ne.nelec).or.
c    &      (nelec.lt.mult-1)) call vberr(' multiplicity does not match 
c    &number of electrons')
c
      else if (label(1:3).eq.'sym') then
c...
c...                   read symmetry designations
c...
         symc = .true.
         isymi = 0
         call input
         call inpi(nbasis)
         call input
         call rdsymt(isymi,nbasis,irconf,group,irrep)
c
      else if (label(1:6).eq.'spin') then
c...
c...                 set the default spinprojection 
c...
         call inpa(label)
         if (label(1:3).eq.'off') then
            projec = .false.
            write(iwr,601)
601         format(/'  *** spinprojection turned off *** ')
         else if (label(1:5).eq.'rumer'.or.label(1:2).eq.'lt') then
            spin =.true.
            write(iwr,602)
602         format(/'  *** produce rumer spinfunctions ***')
         else if (label(1:2).eq.'bd'.or.label(1:2).eq.'yk'.or.
     *            label(1:3).eq.'gen') then
            spin = .false.
            write(iwr,603)
603         format(/'  *** produce branching diagram spinfunctions ***')
         else            
            write(iwr,33) label
            call vberr(' unrecognised spin-scheme')
         end if
c    
      else if (label(1:4).eq.'stru'.or.
     &         label(1:4).eq.'conf') then
         call inpa(label1)
         if ( label(1:4).eq.'stru'.or.
     &       (label(1:4).eq.'conf'.and.label1(1:3).eq.'man')) then
         if (label1(1:3).eq.'man') call inpa(label1)
         onorm = .false.
         if (label1(1:4).eq.'norm') onorm = .true.
c...
c...                   read ready-made structures
c...
         manual = .true.    
         ntdet  = 0
         nconf = 0
         ntdeto = 0
c
51       call input
         call inpa(label)
         call cinput(jrec,jump,-1)
         if (label(1:4).eq.'stru'.or.label(1:3).eq.'end') then
            if (ntdet.eq.0) call caserr('no first structure')
            nconf = nconf + 1
            nspath(nconf) = ntdet - ntdeto
            ntdeto = ntdet
            if (onorm) then
               anorm = 0.0d0
               do i=ntdet-nspath(nconf)+1,ntdet
                  anorm = anorm + coefff(i)*coefff(i)
               end do
               anorm = 1.0d0/dsqrt(anorm)
               do i=ntdet-nspath(nconf)+1,ntdet
                  coefff(i) = coefff(i)*anorm
               end do
            end if
            if (label(1:4).eq.'stru') then 
               call inpa(label1)
               call inpa(label1)
               onorm = .false.
               if (label1(1:4).eq.'norm') onorm = .true.
               call input
            else
               go to 1
            end if
            if (nconf.gt.mconf) call dimerr('configs',mconf)
         end if
         if (label(1:3).ne.'end') then
            ntdet = ntdet + 1
            if (ntdet.gt.mdet) call dimerr('determinants',mdet)
            call inpf(coefff(ntdet))
            call rdista(kdeter(1,ntdet),numb,nelec+2*ncore_i)
c
            if (ncore_i.gt.0) then
               numb = numb - 2*ncore_i
               do i=1,numb
                  kdeter(i,ntdet) = kdeter(i+ncore_i*2,ntdet) - ncore_i
               end do
            end if
            if (numb.ne.nelec) call vberr(' wrong number of electrons in
     & current determinant')
c...      get in correct aaaaaa bbbbb order doubles do not matter for sign !!
c            ndouble = 0
c            do i=1,nelec,2
c               if (kdeter(i,ntdet).eq.kdeter(i+1,ntdet)) then
c                  ndouble = ndouble + 1
c                  if (ndouble.ne.i/2+1) call caserr('doubles not first')
c               end if
c            end do
c            print *,'ndouble: ',ndouble
c...        get in right aaaaabbbb order
c            nalpha = (nelec+mult-1)/2
c            nbeta  = nelec - nalpha
c            if (ndouble.gt.0) then
c               ka = 0
c               kb = nalpha
c               do i=1,nalpha-ndouble
c                  kdeter(i,ntdet+1) = kdeter(ndouble*2+i,ntdet)
c               end do
c               do i=1,nbeta-ndouble
c                  kdeter(i+nalpha,ntdet+1) = 
c     1             kdeter(ndouble*2+nalpha-ndouble+i,ntdet)
c               end do
c               do i=1,ndouble
c                 kdeter(nalpha-ndouble+i,ntdet+1) = kdeter(i*2,ntdet)
c                 kdeter(nalpha+nbeta-ndouble+i,ntdet+1) = 
c     1                 kdeter(i*2,ntdet)
c               end do
c               do i=1,nelec
c                  kdeter(i,ntdet) = kdeter(i,ntdet+1)
c               end do
c...        we move ndouble orbitals through nbeta orbitals
c            if ((ndouble/2)*2.ne.ndouble.and.(nbeta/2)*nbeta.ne.nbeta)
c     1         coefff(ntdet) = -coefff(ntdet)
c            end if
c...
            go to 51
         end if  
c
         else 
c...
c...                 read ready-made configurations
c...
         if (mrsd) call vberr(' no configurations outside the reference 
     &space may be specified')
c
53       call input
         call inpa(label)
         call cinput(jrec,jump,-1)
         if (label(1:3).ne.'end') then
c
c...     conf all directive
c...     generates all possible combinations and occupancies for 'n' orbitals
c...     and 'm' electrons.
c
           if ( label(1:3) .eq. 'all' ) then
             all_nrestr = 0
             all_covalent = 0
             all_ionic = 0
             all_noionic = 0
             all_rumer = .false.
             call inpa(label)
             call inpi(all_norb)
             call inpi(all_nel)
             write(iwr,*)
             write(iwr,'(a,I3,a,I3,a)') 
     &       ' Distributing ',all_nel,' electrons over ',all_norb,
     &       ' orbitals. '
             write(iwr,'(a)') ' Generating all possible combinations. '
531          call inpa(label)
             if ( label(1:5) .eq. 'restr' ) then
               call inpi(all_nrestr)
               do i=1,all_nrestr
                 call inpi(all_restr(i))
               end do
               write(iwr,'(a,30I4)') ' Restricting orbitals: ',
     &         ( all_restr(i), i=1,all_nrestr )
               goto 531
             else if ( label(1:4) .eq. 'allc' ) then
               all_covalent = 1
               write(iwr,'(a,a)') 
     &         ' Generating all possible structures for each covalent ',
     &         'configuration.'
               goto 531
             else if ( label(1:4) .eq. 'alli' ) then
               all_ionic = 1
               write(iwr,'(a,a)') 
     &         ' Generating all possible structures for each ionic ',
     &         'configuration.'
               goto 531
             else if ( label(1:3) .eq. 'noi' ) then
               all_noionic = 1
               write(iwr,'(a,a)')
     &         ' No ionic configuration will be generated. '
               goto 531
             else if ( label(1:3) .eq. 'rum' ) then
               all_rumer = .true.
               write(iwr,'(a)')
     &         ' RUMER diagrams used to generate spin-functions. '
               goto 531
             else
               call cinput(jrec,jump,-1)
             end if
             call allconf_drv(kconf,kscr,nconf,rumer,ispath,nspath,
     &                        nelec,mconf,melec1,maxsp)
             goto 53
           end if
            nconf = nconf + 1
            if (nconf.gt.min(maxcon,mconf)) call dimerr('configurations
     &',min(maxcon,mconf))
            call rdista(kscr,numb,nelec+2*ncore_i)
            call inpa(label)
            if (label(1:3).eq.'kek') then
               call rdista(kscr(numb+1),numkek,nelec+2*ncore_i-numb)
               numb = numb + numkek
            end if
            if (ncore_i.gt.0) then
               numb = numb - 2*ncore_i
               do i=1,numb
                  kscr(i) = kscr(i+ncore_i*2) - ncore_i
               end do
            end if
            if (numb.ne.nelec) call vberr(' wrong number of electrons in
     & current configuration')   
            if (label(1:5).eq.'rumer'.or.label(1:2).eq.'lt') then
               rumer(nconf) = .true.
               call rdist1(ispath(1,nconf),nspath(nconf))
               if (nspath(nconf).gt.maxsp) call dimerr('selected spin pa
     &ths',maxsp)
            else if (label(1:3).eq.'kek') then
               fdist = 0.0d0
               call inpa(label1)
               if (label1(1:3).eq.'dis') call inpf(fdist)
               call kekule(v,kconf,kscr,mconf,melec1,ncore_i,nelec,
     &                     ispath,nspath,numkek,maxsp,nconf,fdist)
            else if (label(1:2).eq.'bd'.or.label(1:2).eq.'yk'.or.
     &               label(1:3).eq.'gen') then
               rumer(nconf) = .false.
               call rdist1(ispath(1,nconf),nspath(nconf))
               if (nspath(nconf).gt.maxsp) call dimerr('selected spin pa
     &ths',maxsp)
            else if (label(1:1).eq.' ') then
               rumer(nconf)  = spin
               nspath(nconf) = 0
            else            
               write(iwr,33) label
               call vberr(' unrecognised spin-scheme')
            end if     
            if (label(1:3).ne.'kek') call reconf(kconf,melec1,nelec,
     &                                           kscr,nconf)
            go to 53
         else
            if (nconf.gt.min(maxcon,mconf))
     &                   call dimerr('configurations',min(maxcon,mconf))
            goto 1
         end if
         endif
c
      else if (label(1:4).eq.'mrsd') then
c...
c...           multi-reference singles and doubles option
c...
         if (nconf.gt.0) call vberr(' input failure in multi-reference s
     &ingles/doubles - configurations defined beside reference')
c...  determine excitation type
         call inpa(label) 
         if (label(1:4).eq.'stru') then
            excstr = .true.
         else if (label(1:4).eq.'conf') then
            excstr = .false.
         else if (label(1:1).eq.' ') then
            call vberr(' input failure in multi-reference singles/double
     &s - no excitation type specified')
         else                                     
            write(iwr,33) label
            call vberr(' input failure in multi-reference singles/double
     &s - unknown excitation type')
         end if
         mrsd   = .true.
         nrefco = 0
         nvirt  = 0               
         next   = .false.
147      if (.not.next) call input 
         call inpa(label)
         next   = .false.
         if (label(1:4).eq.'refe') then  
c...
c...  gather reference configurations    
c...                         
            if (nrefco.gt.0) call vberr(' input failure in multi-referen
     &ce singles/doubles - multiple definition of reference')
            messag = .false.
            nsp0   = 0
            call input
c...  read a reference-configuration on scratch
144         call rdista(kscr,numb,nelec+2*ncore_i)
            if (ncore_i.gt.0) then
               numb = numb - 2*ncore_i
               do i=1,numb
                  kscr(i) = kscr(i+ncore_i*2) - ncore_i
               end do
            end if
            if (numb.ne.nelec) call vberr(' input failure in multi-refer
     &ence singles/doubles - wrong number of electrons in reference conf
     &iguration')
c92            if (nrefco+1.gt.min(maxcon,mconf)) call dimerr ('configurat
c92     &ions',mconf)
            if (nrefco+1.gt.min(maxcon,mconf)) call dimerr ('configurati
     &ons',min(maxcon,mconf))
c...  look for spin-forbidden configurations
            nsocc = 0
            do 146 iorb1=1,nelec
               do 148 iorb2=1,nelec
                  if (iorb2.eq.iorb1) go to 148
                  if (kscr(iorb1).eq.kscr(iorb2)) go to 146
148            continue
               nsocc = nsocc + 1
146         continue
            if (nsocc.lt.mult-1) then
               if (excstr) then
                  call vberr(' input failure in multi-reference singles/
     &doubles - a reference configuration is spin forbidden')
               else 
                  write(iwr,355)
355               format(/10x,'warning!!!!!'//,
     &                 '-- a reference configuration is spin forbidden')
               end if
            end if
c...  pick up requests for structure selection/alternative spin
c...  projection. forbid this for excitations on configuration basis
c...  (configurations with and without structure selection are put at
c...  different locations on kconf)
            call inpa(label)  
            call inpi(idum)
            if (idum.eq.0) then
c...  no structure selection: configurations without structure selection
c...  (nspath=0) must appear first on kconf, so make room for current 
c...  configuration
               if (nrefco.gt.nsp0) then
                  if (.not.messag) then
                     write(iwr,435) 
435                  format(/10x,'message'//,
     &                     '-- configurations with structure ',  
     &                       'selection are put behind others')
                     messag = .true.    
                  end if
                  call copi( kconf(1,nsp0+1), kconf(1,nsp0+2),
     &                      (nrefco-nsp0)*melec1,-1)
                  call copl( rumer(nsp0+1), rumer(nsp0+2),
     &                      nrefco-nsp0,-1)
                  call copi(nspath(nsp0+1),nspath(nsp0+2),
     &                      nrefco-nsp0,-1)
                  call copi(ispath(1,nsp0+1),ispath(1,nsp0+2),
     &                      (nrefco-nsp0)*maxsp ,-1)
               end if
               jconf = nsp0 + 1
               nsp0  = nsp0 + 1  
            else
               if (.not.excstr) call vberr(' input failure in multi-refe
     &rence singles/doubles - illegal selection of structures from refer
     &ence configuration') 
c...  put current configuration last
               jconf = nrefco + 1  
               call cinput(jrec,jump,-1)
            end if
c...  store current configuration on column jconf of kconf
c...  put doubles last, order them and store their number at row nelec+1
            call reconf(kconf,melec1,nelec,kscr,jconf)
c...  change spin defaults if requested
            if (label(1:5).eq.'rumer'.or.label(1:2).eq.'lt') then
               rumer(jconf) = .true.
               call rdist1(ispath(1,jconf),nspath(jconf))
               if (nspath(jconf).gt.maxsp) call dimerr('selected spin pa
     &ths',maxsp)
            else if (label(1:2).eq.'bd'.or.label(1:2).eq.'yk'.or.
     &               label(1:3).eq.'gen') then
               rumer(jconf) = .false.
               call rdist1(ispath(1,jconf),nspath(jconf))
               if (nspath(jconf).gt.maxsp) call dimerr('selected spin pa
     &ths',maxsp)
            else if (label(1:1).eq.' ') then
c...  no structure selection: order singles as well (their relative
c...  positions have no significance)
               rumer(jconf)  = spin
               nspath(jconf) = 0
           call ordernumbs(kconf(1,jconf),nelec-2*kconf(nelec+1,jconf))
            else            
               write(iwr,33) label
               call vberr(' input failure in multi-reference singles/dou
     &bles - unrecognised spin-scheme')
            end if
c...  current configuration-definition complete
            nrefco = nrefco + 1    
            call input         
            call inpa(label)  
            call cinput(jrec,jump,-1)
            if (locatc(figure,10,label(1:1)).gt.0) then 
               go to 144 
            else            
               next = .true.
            end if
c...  collect internal orbitals and order them
            call colint(kconf,nelec+1,nrefco,korb,ninner)
            if (korb(ninner).gt.253) call vberr(' input failure in multi
     &-reference singles/doubles -  highest orbital number in reference 
     &space exceeds 253')
c...  add symbolic externals
            if (n8_16.eq.0) call caserr('choose 8/6 bits packing')
            korb(ninner+1) = 2**n8_16-2
            korb(ninner+2) = 2**n8_16-1
            norb = ninner + 2              
c...  store configuration identity and remove redundant configurations
            do 440 iconf=1,nrefco
               call storid(kconf(1,iconf),iconf,idconf,korb,norb)   
440         continue 
            nrejco = 0
            jfirst = 1
450         continue
c...  find configurations equal to jfirst 
            jlast  = jfirst
460         if (jlast+1.le.nrefco-nrejco) then
               if (idconf(3,jlast+1).eq.idconf(3,jfirst)) then
                  if (idconf(2,jlast+1).eq.idconf(2,jfirst)) then
                     if (idconf(1,jlast+1).eq.idconf(1,jfirst)) then
c...  configuration jlast+1 equals jfirst
                        jlast = jlast + 1
                        go to 460
                     end if
                  end if
               end if
            end if    
            nnrejc = 0
            if (jlast.gt.jfirst) then  
               if (nspath(idconf(jfirst,4)).eq.0) then  
c...  since no structure selection is applied to first, all other are 
c..   redundant
                  do 470 iconf=jfirst+1,jlast     
                     write(iwr,465) idconf(4,iconf)
465                  format(/10x,'warning!!!!!'//,
     &                  '-- configuration', i3, ' is redundant')    
                     kscr  (nnrejc+1) = iconf             
                     nnrejc = nnrejc + 1
                     krejco(nrejco+1) = idconf(4,iconf)
                     nrejco = nrejco + 1    
470               continue  
                  call skipi(idconf,4,nrefco-(nrejco-nnrejc),kscr,
     &                       nnrejc,idconf)
               end if   
            end if
            jfirst = jlast - nnrejc + 1
            if (jfirst.lt.nrefco-nrejco)  go to 450
            if (nrejco.gt.0) then          
c...  shift configuration definitions
               call skipi(kconf ,melec1,nrefco,krejco,nrejco,kconf)
               call skipl(rumer ,1     ,nrefco,krejco,nrejco,rumer)
               call skipi(nspath,1     ,nrefco,krejco,nrejco,nspath)
               call skipi(ispath,maxsp ,nrefco,krejco,nrejco,ispath)
c...  correct explicit position configuration on idconf
               do 490 iconf=1,nrefco
                  jshift = 0
510               if (idconf(4,iconf).gt.krejco(jshift+1)) then
                     jshift = jshift + 1
                     go to 510
                  end if
                  idconf(4,iconf) = idconf(4,iconf) - jshift
490            continue
               nrefco = nrefco - nrejco   
            end if  
            nconf  = nrefco    
            go to 147
         else if (label(1:4).eq.'virt') then
c...     read virtual orbitals                         
            if (nvirt.gt.0) call vberr(' input failure in multi-referenc
     &e singles/doubles - virtuals already defined')
            call input
            call rdista(kvirt,nvirt,mxorbvb)
            call ordernumbs(kvirt,nvirt)                      
            if (n8_16.eq.0) call caserr('choose 8/6 bits packing')
            if (kvirt(nvirt).gt.2**n8_16-1) 
     & call vberr('  input failure in multi-reference singles/doubles - 
     &highest orbital number in virtual space exceeds maximum')
            if (nvirt.le.1)  call vberr(' input failure in multi-referen
     &ce singles/doubles - not sufficient virtuals specified')     
            go to 147
c                                           
         else if (label(1:3).eq.'end') then
            if (nrefco.eq.0)        call vberr(' input failure in multi-
     &reference singles/doubles - no reference specified')     
            if (nvirt.eq.0)         call vberr(' input failure in multi-
     &reference singles/doubles - no virtuals specified')
         else 
            write(iwr,33) label
            call vberr(' input failure in multi-reference singles/double
     &s - unrecognised directive')
         end if
c
c...  =aivb=
c...  choosing atomic states, automatic generation of determinants
c...  08-09/2007, 04-11/2008, 03-11/2009, 05/2010 09-10/2010 - Marcin
c
      else if ( label(1:4) .eq. 'aivb' ) then
c
c...  read atomic state definitions
c
      write(iwr,7110)
7110  format(/8x,' %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ',
     1       /8x,' #                                                # ',
     2       /8x,' #          Atoms in Valence Bond (AiVB)          # ',
     3       /8x,' #               version March 2012               # ',
     4       /8x,' #                                                # ',
     5       /8x,' #        email: Marcin.Zielinski@sara.nl         # ',
     6       /8x,' #                                                # ',
     7       /8x,' %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ')
         manual = .true.
         nread = 0
         ntdet  = 0
         nconf = 0
         nconfp = 0
         ntdetp = 0
         ntdeto = 0
         aivb_reduce = .false.
         aivb_diff = .false.
         aivb_set = .true.
         aivb_moperm = .false.
         aivb_optorb = .false.

         do iijj=1,14
           aimo_debugpr(iijj) = 0
         end do
c
         write(iwr,'(a)') ' '
71       call input
         call inpa(label)
         if ( label(1:4) .eq. 'redu' ) then 
           aivb_reduce = .true.
           write(iwr,'(a,a)') ' => AiVB removing redundant ',
     &     'configurations requested <='
           call vberr('AiVB reduce does not work properly at the moment
     &- do not use it.') 
           call inpf(fthr)
           aivb_threshold = fthr
           if ( aivb_threshold .eq. 0.0d0 ) then 
             call cinput(jrec,jump,-1)
             aivb_threshold = 10E-08
           end if
           write(iwr,'(a,E10.4)') ' => AiVB Threshold: ',
     &                                aivb_threshold
           goto 71
         else if ( label(1:4) .eq. 'diff' ) then
           aivb_diff = .true.
           write(iwr,'(a,a)') ' => AiVB different atomic orbitals ', 
     &     'for different atomic states requested <='
           goto 71
         else if ( label(1:4) .eq. 'debu' ) then
           write(iwr,'(a)') ' => AiVB debug printouts requested '
           do iijj=1,14
             call inpi(aimo_debugpr(iijj))
           end do
           goto 71
         else if ( label(1:4) .eq. 'orbo' ) then
           aivb_optorb = .true.
           write(iwr,'(a)') ' => AiVB orbital optimization ON <= '
           goto 71
         else if ( label(1:4) .eq. 'mope' ) then
           aivb_moperm = .true.
           call rdist1(aivb_mo2perm,aivb_nmo2p)
           write(iwr,'(a,10(I3,a,I3,a))') 
     &     ' => AiVB permuting orbitals: ',
     &     (aivb_mo2perm(iijj),' with ',aivb_mo2perm(iijj+1),' - ', 
     &     iijj=1,aivb_nmo2p,2)
           goto 71
         else if ( label(1:4) .eq. 'conf' ) then
           goto 709
         else if ( label(1:3) .eq. 'end' ) then
           if ( nconf .eq. 0 ) then
             call vberr('no configurations specified within AIMO CONF!')
           end if
           goto 708
         end if

709      continue

         nread = nread + 1 
         if ( aimo_debugpr(1) .gt. 0 ) then
           write(iwr,*) ' ====> READ beginning values '
           write(*,'(a,I3,a,I3)') ' nread   = ',nread,
     1     '  ntdet = ',ntdet
           write(*,'(a,I3,a,I3)') ' ncore_i = ',ncore_i,
     2     '  nelec = ',nelec
           write(iwr,'(a,I10)') ' mdet = ',mdet
         end if
         if (ntdet.gt.mdet) call dimerr('determinants',mdet)
c
c...  read input definitions
c
         call read_aivb_input(nread)
         aimo_confnr(nread) = 0
c
c...  generate atomic state configurations and determinants
c
         call generate_aivb_drv(coefff,kdeter,nread,ntdet,nconf,numb,
     &        nelec+2*ncore_i,nelec,ncore_i,mdet,kscr)
c
         if ( aimo_debugpr(1) .gt. 0 ) then
           write(*,'(a,5I4)') ' 1 ncore_i,nelec,numb,ntdet,nconf: ',
     1                         ncore_i,nelec,numb,ntdet,nconf
         end if

         numb = numb - 2*ncore_i

         if ( aimo_debugpr(1) .gt. 0 ) then
           write(*,'(a,5I4)') ' 2 ncore_i,nelec,numb,ntdet,nconf: ',
     1                         ncore_i,nelec,numb,ntdet,nconf
         end if

         if ( aimo_debugpr(1) .gt. 0 )
     1     write(iwr,6203) 
6203     format(70('-'))

         if (numb.ne.nelec) call vberr(' wrong number of electrons in
     & current determinant')
         ntdeto = ntdet
c...
         go to 71

708      continue

         if (aivb_moperm) then
           do k=1,aivb_nmo2p,2
             do i=1,ntdet
               do j=1,numb
                 if ( kdeter(j,i) .eq. aivb_mo2perm(k) ) then
                   if ( aimo_debugpr(1) .gt. 1 ) then
                     write(iwr,'(a,I3,a,a,I3,a,I3)') 
     &               ' det nr. ',i,': ',
     &               ' replacing orb. ',kdeter(j,i),
     &               ' with ',aivb_mo2perm(k+1)
                   end if
                   kdeter(j,i) = aivb_mo2perm(k+1)
c                   coefff(i) = coefff(i)*(-1.0d0)
                 else if ( kdeter(j,i) .eq. aivb_mo2perm(k+1) ) then
                   if ( aimo_debugpr(1) .gt. 1 ) then
                     write(iwr,'(a,I3,a,a,I3,a,I3)') 
     &               ' det nr. ',i,': ',
     &               ' replacing orb. ',kdeter(j,i),
     &               ' with ',aivb_mo2perm(k)
                   end if
                   kdeter(j,i) = aivb_mo2perm(k)
c                   coefff(i) = coefff(i)*(-1.0d0)
                 end if
               end do
             end do
           end do
         end if

         if ( aimo_debugpr(1) .gt. 1 ) then
           do j=1,ntdet
             write(*,'(1x,F7.4,a,30I4)') coefff(j),' : ',
     1            ( kdeter(i,j), i=1,numb)
           end do

           write(iwr,'(a)') 'determinants for STRUC: '
           do j=1,ntdet
             write(*,'(1x,F7.4,a,I3,a,30I4)') 
     2            coefff(j),' core ',ncore_i,' ',
     1            ( kdeter(i,j)+ncore_i, i=1,numb)
           end do

         end if

         if (ntdet.eq.0) call caserr('no first structure')
c...  loop over the number of input reads 'nread'
         do jj=1,nread
c...  loop over the number of configurations per 'nread'
           do ii=1,aimo_confpr(jj)
             nspath(ii+nconfp) = aimo_detstruc(jj,2,ii)
             ntdetp = ntdetp + aimo_detstruc(jj,2,ii)
             if ( aimo_debugpr(1) .gt. 0 ) then
               write(iwr,'(1x,a,a,5I5)') 'AiVB: ',
     &         'nspath,nconfp,ntdetp',
     &         nspath(ii+nconfp),nconfp,ntdetp
             end if
           end do
           nconfp = nconfp + aimo_confpr(jj)
         end do
         if (nconf.gt.mconf) call dimerr('configs',mconf)
         if ( aimo_debugpr(1) .gt. 0 ) then
           write(iwr,'(a,2I5)') ' nconf,mconf = ',nconf,mconf
           print *,'label after: ',label(1:4)
         end if
      else if (label(1:4).eq.'crit') then
         call inpf(xxx)
         call inpi(i) 
         xxx = xxx * 10**(i*1.0d0)
         cridep = xxx
         write(iwr,6204) cridep
6204     format(/,' dependency criterion ',e12.5)
      else
         if (label(1:3).ne.'end') then 
            write(iwr,33) label
            call vberr(' unrecognised directive')
         end if
      end if
c
      if (label(1:4).ne.'end') go to 1
c...
c...  ************************ end cycle *******************************
c...
      write(iwr,6999) n8_16,2**(n8_16)-1,mxorbvb
6999  format(/8x,' $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$',
     1       /8x,' $$           ',i6,'-bits packing in VB           $$',
     2       /8x,' $$           maximum # orbitals :',i6,'          $$',
     3       /8x,' $$     max. # orbitals in turtleparam :',i6,'    $$',
     4       /8x,' $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$')
c
      if (n8_16.eq.8) then
         n8_16 = 8
         n16_32 = n8_16*2
         n340 = 340
c..         kup1 = 341         n340+1
c..         kup2 = 341+85      n340+1+n340/(32/n8_16)
      else if (n8_16.eq.16) then
         n8_16 = 16
         n16_32 = n8_16*2
         n340 = 254
c..         kup1 = 255            n340+1
c..         kup2 = 255 + 127      n340+1+n340/(32/n8_16)
      else
         call caserr('unexpected n8_16')
      end if
      if (n340.gt.340) call caserr(' vbtran overflow ')
c
      if (mult.eq.0) then
c...
c...  set default multiplicity
c...
         if ((nelec/2)*2.eq.nelec) then
            mult = 1
         else
            mult = 2
         end if
      end if
c
      if (manual) return  
c
      if (.not.mrsd) then
c...
c...  reject non-unique configurations
c...
         call redcon(kconf,kscr,krejco,melec1,nconf,nelec,nrejco,nspath)
         do 480 irejco=1,nrejco
            write(iwr,482) krejco(irejco)
482         format(/10x,'warning!!!!!'//,
     &      '-- configuration', i3, ' is redundant')
480      continue
         nconf = nconf - nrejco
      end if
c
      if (allout) then
c...
c...  write configuration-info to output
c...
            write(iwr,25)
25          format(/,' ',14x,'**** all unique configurations **** ')
            object = '  configuration'
         call writ2d(kconf,melec1,object,nelec,1,nconf)
      endif
c
      if (symc) then
c...
c...  reject symmetry-forbidden configurations
c...
         nrejco = 0
         do 250 iconf=1,nconf
            isym = isconf(kconf(1,iconf),nelec)
            if (isym.ne.irconf) then
               krejco(nrejco+1) = iconf
               nrejco = nrejco + 1
            endif
250      continue
         if (nrejco.gt.0) call skipi(kconf,melec1,nconf,krejco,nrejco,kc
     &onf)
         nconf = nconf - nrejco
      endif
c...
c...  calculate maximum number of open shells and reject spin-forbidden 
c...  configurations (this has been done already for reference 
c...  configurations)
c...
      msocc = 0
      nrejsp = 0
      if (mrsd) then
         do 360 iconf=1,nconf
            msocc = max(nelec-2*kconf(melec1,iconf),msocc)
360      continue
      else
         do 350 iconf=1,nconf
            msocc = max(nelec-2*kconf(melec1,iconf),msocc)
            if (nelec-2*kconf(melec1,iconf).lt.mult-1) then
               krejco(nrejsp+1) = iconf
               nrejsp = nrejsp + 1
            end if
350      continue
         if (nrejsp.gt.0) then
            call skipi(kconf ,melec1,nconf,krejco,nrejsp,kconf)
            call skipl(rumer ,     1,nconf,krejco,nrejsp,rumer)
            call skipi(nspath,     1,nconf,krejco,nrejsp,nspath)
            call skipi(ispath,maxsp ,nconf,krejco,nrejsp,ispath)
            nconf = nconf - nrejsp
         end if
      end if
c
      if (nconf.le.0) call vberr(' no accurate configurations specified'
     &)
c...
c...  end of input/initialize vital variables
c...
c.ab
c.ab    set highest spin path as defined by input
c.ab
      msp = 0
      nsocchi = 0
      do 502 iconf=1,nconf
         ndocci = 0
         do ielec=1,nelec-1
          if ( kconf(ielec,iconf) .eq. kconf(ielec+1,iconf) ) then
            ndocci = ndocci + 1
          end if
         end do
         nsocci = nelec - ndocci*2
         if ( nsocci .gt. nsocchi ) then
           nsocchi = nsocci
           if ( nsocci .lt. 3 ) then
             msp = 1
           else
             spinmax = (mult-1.0d0)/2.0d0
             n2s = nsocci/2.00d0 - spinmax
             new1 = newtont(nsocci,n2s)
             n2s1 = n2s-1
             new2 = newtont(nsocci,n2s1)
             mspmax = new1 - new2
             do 501 isp=1,maxsp
                msp = max(ispath(isp,iconf),mspmax)
501          continue
           end if
         end if
502   continue
c.ab
      call initsy(kconf,melec1,nconf,nspath,msp)
c
      write(iwr,135)
      if (symc) then
         write(iwr,145) group
         write(iwr,155) irrep
         if (nrejco.gt.0) write(iwr,165) nrejco
      endif
      write(iwr,175) mult
      write(iwr,185) nelec
      if (nrejsp.gt.0) write(iwr,195) nrejsp
      write(iwr,205) nconf
      write(iwr,215)
135   format(/1x,45(1h*),/,t11,'system definitions')
145   format(              t10,'point-group              : ',a)
155   format(              t10,'irrep. of wave function  : ',a)
165   format(              t10,'rejected on base of sym. : ',i6,
     &                                               ' configurations')
175   format(              t10,'multiplicity             : ',i6)
185   format(              t10,'number of electrons      : ',i6)
195   format(              t10,'rejected on base of spin : ',i6,
     &                                               ' configurations')
205   format(              t10,'number of configurations : ',i6)
215   format( 1x,45(1h*))
c
      object = 'configuration'
      call writ2d(kconf,melec1,object,nelec,1,nconf)
c...
c...  store configuration info on ed7
c...
      k7conf = kscra7vb('k7conf',nconf*(nelec+1),'i','w')
      call wrt3is(kconf,(nelec+1)*nconf,num8)
c
      return
      end
      subroutine crestr(q,lword,v,qq)
c...
c...   main program for symbolic vb calculation specification
c...   j.h. langenberg, j. verbeek, j.h. van lenthe
c...   theoretical chemistry group
c...   utrecht 1988/1989/1990
c...
c
c     'crestr' can perform all possible single and double excitations of
c     a given reference-function to a given virtual space.
c     it can also handle hand-given configurations, or a combination of
c     these with generated configurations.
c     with each configuration, a number of structures is made. each
c     structure represents one of the spin-states with multiplicity
c     '2s+1', that can be accomplished with the actual number of single
c     occupied orbitals. the generation of structures is based on the
c     'leading terms'-scheme.
c     'crestr' is able to skip configurations on base of symmetry.
c
c
      implicit REAL (a-h,o-z), integer (i-n)
c
INCLUDE(common/turtleparam)
INCLUDE(common/c8_16vb)
c
      common  /indep/ time0,nelec,mult
      common   /hans/ msocc,mstrco,mdetco,mdetst,mperm,mstruc
INCLUDE(common/logics)
      common /size/   nrefco,nsinin,nsinex,ndouin,ndoumi,ndoue1,ndoue2
INCLUDE(../m4/common/iofile)
c
      character*15 object
      character* 8 test
c
      common/crestn/ struct
INCLUDE(../m4/common/work)
INCLUDE(../m4/common/timeperiods)
INCLUDE(common/splice)
INCLUDE(common/first_vb)
INCLUDE(common/tractlt)
      character*44 struct,dummy
      dimension q(*),v(*),qq(*)
c
c     data maxsp/14/
c
c     data mstrco,mdetco,mdetst,mperm calculated in initsy
c

      maxsp = maxspvb
      time0 = cpulft(1)
cfd     print *,'TP_VB_STRUC',TP_VB_STRUC
      call start_time_period(TP_VB_STRUC)
_IF(atmol)
      dummy = ''
_IF(parallel)
      write(dummy,'(A1I4.4)') '.',ipg_nodeid()
_ENDIF
      struct = 'structures'//dummy(1:6)
      write(iwr,*) 'structure file',struct
      call strtrm(struct,nn)
      open(15,form='unformatted',file=struct(1:nn))
_ELSE
      struct = 'on ed7'
_ENDIF
      write(iwr,12) struct
12    format(' && structures-file ',a,' &&',/)
c
      ncore = 0
      ncore_i = 0
      ofirst = .true.
      opr_first = .false.
C
C*** TO BE SUPPLIED BY GAMESS
C
      call gam_to_vb(nelec,mult,nbasis)
      if (nbasis.le.255) then
         n8_16 = 8
      else
         n8_16 = 16
      end if
c...  read number of electrons
c      call input
c10    call inpa(test)
c      if (test.eq.'nelec') then
c         call inpi(nelec)
c         go to 100
c      else if (test.eq.' ') then
c         call input
c      else
c        call vberr(' unrecognised post-predirective for symbolic')
c      end if
c      go to 10

c...  Check if the number of electrons passed by GAMESS should be overruled...
20    call input
      call inpa(test)
c...
c...  if so, then change it and read the next line
      if (test.eq.'nelec') then
c...
c...  override the number of electrons supplied by GAMESS
         call inpi(nelec)
         write(iwr,*) "Overriding number of electrons with ",nelec," 
     +            electrons"
         goto 20
      else if (test.eq.'mult') then
c...
c...  override the multiplicity supplied by GAMESS
         call inpi(mult)
         write(iwr,*) "Overriding the multiplicity with",mult
         goto 20
      else if (test.eq.'core') then
c...     read in number of core orbitals and adjust nelec
         call inpi(ncore_i)
         write(iwr,*) ' # orbitals treated as core : ',ncore_i
         nelec = nelec - 2*ncore_i
         ncore = ncore_i
         do i=1,ncore
           mapcie(i) = i
         end do
         goto 20
      else if (test(1:4).eq.'255-'.or.test(1:4).eq.'<255') then
         n8_16 = 8
         goto 20
      else if (test(1:4).eq.'255+'.or.test(1:4).eq.'>255') then
         n8_16 = 16
         goto 20
      else 
c...  If not so, back one position on current input-line.
         jrec=jrec-1
      end if 
      melec1=nelec+1
100   continue  
c...
c...   read input  / set up core-partitioning for that
c...
      mconf = (lword-mxorbvb-mxorbvb)/
     &        (melec1+1+1+maxsp+4+1+melec1+maxsp*nelec*melec1)
      write(iwr,600) nelec,mconf,lword
600   format(//11x,20(1h*),/,11x,'  vb - symbolic ',/,11x,20(1h*),////,
     &       5x,' parameters  for input.... ',/,5x,' nelec  ',
     &       i6,4x,' mconf  ',i8,4x,' lword  ',i10)

      kconf  = 1 
      nspath = kconf  + mconf*melec1
      ispath = nspath + mconf
      korb   = ispath + mconf*maxsp
      kvirt  = korb   + mxorbvb  
      idconf = kvirt  + mxorbvb
      krejco = idconf + mconf*4
      kkscr  = krejco + mconf


c...  kend   = kkscr  + mconf*melec1 
c.ab  write(iwr,*) 'kconf,nspath,ispath,korb,kvirt,idconf'                
c.ab  write(iwr,*)  kconf,nspath,ispath,korb,kvirt,idconf
c...  partitioning in case of 'structures manual'      
      mdet = lword/(1+nelec)              
c                              
      kcofff = 1
      kdeter = kcofff + mdet
c...  kend   = mdet*nelec
c
      call cresti(q(kconf),melec1,q(nspath),q(ispath),maxsp,
     &            q(korb),q(kvirt),q(idconf),q(krejco),q(kkscr),mconf,
     &            nconf,nsp0,ninner,nvirt,q(kcofff),q(kdeter),
     &            mdet,ntdet,v,qq)
c...        
c...   set up core partitioning for symbolic now we know the problem
c...         
cccccc
cccccc alleen nog eventueel aangegeven verbeteringen overigen      
cccccc integers helft van de ruimte geven: (...+1)/2 <-- nav !!!!!!!!!!
cccccc       
      mbeta = (msocc - (mult-1))/2
      mpair = mbeta   
      mspine = mstrco*msocc + mstrco*msocc + mpair+1 +   
     &         mpair*mperm*(mpair+1)*mstrco + mdetco + mpair*mperm*mpair 
c   
cremco      if (.not.projec) then      
cremcocccccc
cremcocccccc  spin off is infantile programmed: a spin projection is performed
cremcocccccc  but completetely neglegted --> do it all at once in crestd
cremcocccccc                         call deter to construct det's? 
cremcocccccc  spin selection now neglected: could be used for det.-selection     
cremcocccccc
cremcoc...  real arrays
cremco         kcofff = 1
cremco         kpacde = kcofff + nconf*mdetco     
cremco         nnn    = (mdetco*nelec-1)/(64/n8_16) + 1
cremco         kcoeff = kpacde + nconf*nnn     
cremcoc...  integer arrays
cremco         kk     = kconf
cremco         kconf  = kcoeff + mdetco*mdetst   
cremco         call copi(q(kk),q(kconf),nconf*melec1,-1)
cremco         kstruc = kconf  + nconf*melec1  
cremco         kdeter = kstruc + mdetco*mdetst  
cremco         klogic = kdeter + mdetco*nelec
cremco         knstru = klogic + mstrco*mdetco
cremco         kntdet = knstru + nconf
cremco         kdetps = kntdet + nconf
cremco         idetps = kdetps + nconf*mdetco   
cremco         kkscr  = idetps + nconf*mdetco      
cremcoc
cremco         kend   = kkscr  + mspine  
cremco         mscr   = lword  - kkscr        
cremco      else if (manual) then 
      if (manual) then 
c...  real arrays
         mconf =  nconf
         maxsp =  1
c
         kcofff = 1       
         kpacde = kcofff + ntdet     
c....     take ample
         nnn    = ntdet*((nelec-1)/(64/n8_16)+1)    
c...  integer arrays
         kk     = kdeter    
         kdeter = kpacde + nnn    
         if (kdeter.le.kk) then    
            incr = +1
         else
            incr = -1
         end if
         call copi(q(kk),q(kdeter),ntdet*nelec,incr)

         knstru = kdeter + ntdet*nelec
         kntdet = knstru + nconf
         kdetps = kntdet + nconf
         kbonds = kdetps + nconf
         idetps = kbonds + nconf*maxsp*nelec
         kend   = idetps + ntdet
c
      else if (mrsd) then   
ccccc origineel in buffer  
c...
c...  we know nconf now / get input arrays neatly alligned for safety                
c...
         nnspat = nspath    
         iispat = ispath
         kkvirt = kvirt
         kkorb  = korb
         iidcon = idconf        
c
c...     kconf  = 1
         nspath = kconf  + nconf*melec1  
         ispath = nspath + nconf              
         kvirt  = ispath + nconf*maxsp       
         korb   = kvirt  + nvirt          
         idconf = korb   + ninner+2        
c
         call copi(q(iidcon),q(idconf),nconf*4    ,-1)
         call copi(q(kkorb) ,q(korb)  ,ninner+2   ,-1)   
         call copi(q(kkvirt),q(kvirt) ,nvirt      ,-1)
         call copi(q(iispat),q(ispath),nconf*maxsp,-1)   
         call copi(q(nnspat),q(nspath),nconf      ,-1)
c 
         kkconf = kconf                               
         nnspat = nspath    
         iispat = ispath
         kkvirt = kvirt
         kkorb  = korb
         iidcon = idconf        
c
         nrefco = nconf
c...  
c...                      r e a l - scratch  
c...  
c...  mrsdd:        2*mstrco*mdetco                 (projc-->schmic) 
c...  singc/doubc:  2*mstrco*mdetco                 (projc-->schmic)      
c...  sings/doubs:    mstrco*mdetco + mstrco*mdetco (excits-->lindep)   
c...  remred:       3*mstrco*mdetco                 (lindep)     
c... 
         nscrre = 3*mstrco*mdetco                 
c...
c...                      i n t e g e r - scratch
c...
c...  mrsdd:        nconf +                            (mrsdd/projc) 
c...                mspine                     (projc-->spinef)
c...  singc/doubc:  nelec+1 + 
c...                mspine                     (projc-->spinef)
c...  sings/doubs:  nelec+1 + mdetco +                 (excits)
c...                mdetco*(mdetco+2) +                (excits)
c...                max(mdetco*msocc/2,mdetco+mdetco, (excits)
c...                     msocc+mstrco+mstrco*mdetco+
c...              &      mdetco*nelec+mdetco*nelec+
c...              &      mdetco+msocc/2)               (excits-->lindep)
c...  remred:       mconf + 
c...                mconf +                            (lindep)
c...                max(nelec+nelec,                  (prepli)
c...                     msocc+mstrco+mstrco*mdetco
c...               &     mdetco*nelec+mdetco*nelec+
c...               &     mdetco+msocc/2)               (lindep)
c... 
         nremre = msocc+mstrco+mstrco*mdetco+mdetco*nelec+mdetco*nelec+
     &            mdetco+msocc/2 
c
         nscrin =             nconf + mspine
         nscrin = max(nscrin,nelec+1+mspine)
         nscrin = max(nscrin,nelec+1+mdetco+mdetco*(mdetco+2)+
     &                        max(mdetco*msocc/2,mdetco+mdetco,nremre))
         nscrin = max(nscrin,max(nelec+nelec,nremre))
c
         nnn    = (mdetco*nelec-1)/(64/n8_16) + 1  
c 
         if (allrum) then       
            mconf = (lword - (mstrco*mdetst+mdetco+nscrre+nconf*melec1+
     &                        nconf+nconf*maxsp+nvirt+ninner+2+
     &                        mstrco*mdetst+mdetco*nelec+nscrin)) /
     &              (mstrco*mdetst+nnn+4+1+1+mdetst+mstrco*mdetst+2)
         else                                                              
            mconf = (lword - (mstrco*mdetst+mdetco+nscrre+nconf*melec1+
     &                        nconf+nconf*maxsp+nvirt+ninner+2+
     &                        mstrco*mdetst+mdetco*nelec+nscrin)) /
     &              (mstrco*mdetco+nnn+4+1+1+mdetst+mstrco*mdetco+2)         
         end if           
         if (mconf.le.nconf) call vberr('insufficient core at main parti
     &tioning in crestr')
c...  real arrays
         kcofff = 1
         if (allrum) then  
            kpacde = kcofff + mconf*mstrco*mdetst   
         else
            kpacde = kcofff + mconf*mstrco*mdetco
         end if                          
         kcoeff = kpacde + mconf*nnn      
         kparit = kcoeff + mstrco*mdetst
         kscr   = kparit + mdetco 
c...  integer arrays        
         kconf  = kscr   + nscrre      
         nspath = kconf  + nconf*melec1  
         ispath = nspath + nconf              
         kvirt  = ispath + nconf*maxsp       
         korb   = kvirt  + nvirt          
         idconf = korb   + ninner+2        
         kstruc = idconf + mconf*4
         kdeter = kstruc + mstrco*mdetst
         klogic = kdeter + mdetco*nelec
         knstru = klogic + mstrco*mdetco
         kntdet = knstru + mconf
         kdetps = kntdet + mconf
         idetps = kdetps + mconf*mdetst  
         if (allrum) then  
            kkscr  = idetps + mconf*mstrco*mdetst   
         else                     
            kkscr  = idetps + mconf*mstrco*mdetco   
         end if       
c          
         kend   = kkscr + 2*mconf + nscrin   
         mscr   = lword - (kkscr-1)  
      write(iwr,*) 'kcofff,kpacde,kcoeff,kparit,kscr,kconf,nspath,ispath,
     &kvirt,korb,idconf,kstruc,kdeter,knstru,kntdet,kdetps,idetps,kkscr'    
      write(iwr,*)kcofff,kpacde,kcoeff,kparit,kscr,kconf,nspath,ispath,
     &kvirt,korb,idconf,kstruc,kdeter,knstru,kntdet,kdetps,idetps,kkscr 
c...
c...  now put input-arrays at new adresses
c... 
         call copi(q(iidcon),q(idconf),nconf*4     ,-1)
         call copi(q(kkorb) ,q(korb)  ,ninner+2    ,-1)   
         call copi(q(kkvirt),q(kvirt) ,nvirt       ,-1)
         call copi(q(iispat),q(ispath),nconf*maxsp ,-1)   
         call copi(q(nnspat),q(nspath),nconf       ,-1)   
         call copi(q(kkconf),q(kconf) ,nconf*melec1,-1)
      else  
c...
c...  we know nconf now / get input arrays neatly alligned for safety               
c...
         nnspat = nspath    
         iispat = ispath
c
c...     kconf  = 1
         nspath = kconf  + nconf*melec1  
         ispath = nspath + nconf              
c
         call copi(q(iispat),q(ispath),nconf*maxsp,-1)   
         call copi(q(nnspat),q(nspath),nconf      ,-1)
c
         kkconf = kconf
         nnspat = nspath
         iispat = ispath
c...  real arrays 
         kcofff = 1
         if (allrum) then
            kpacde = kcofff + mstruc*mdetst   
         else
            kpacde = kcofff + mstruc*mdetco
         end if                          
         nnn    = (mdetco*nelec-1)/(64/n8_16) + 1   
         kcoeff = kpacde + nconf*nnn        
         kscr   = kcoeff + mstrco*mdetst
c...  integer arrays    
         kconf  = kscr   + 2*mstrco*mdetco      
         nspath = kconf  + nconf*melec1  
         ispath = nspath + nconf              
         kstruc = ispath + nconf*maxsp
         kdeter = kstruc + mstrco*mdetst
         klogic = kdeter + mdetco*nelec
         knstru = klogic + mstrco*mdetco
         kntdet = knstru + nconf
         kdetps = kntdet + nconf
         kbonds = kdetps + mstruc
         idetps = kbonds + nconf*maxsp*nelec
         if (allrum) then
            kkscr  = idetps + mstruc*mdetst
         else
            kkscr  = idetps + mstruc*mdetco
         end if                  
c
         kend   = kkscr  + mspine
         mscr   = lword  - kkscr                           
         if (mscr.lt.0) then
            write(iwr,610) mscr,lword,kkscr
610         format(/' *** # words available for scratch ',i12,
     1             /' *** # words really available      ',i12,
     2             /' *** start of scratch    at        ',i12)
            call vberr(' enlarge core or reduce problem size')
         end if
c...
c...  now put input-arrays at new adresses
c... 
         call copi(q(iispat),q(ispath),nconf*maxsp ,-1)   
         call copi(q(nnspat),q(nspath),nconf       ,-1)
         call copi(q(kkconf),q(kconf) ,nconf*melec1,-1)
      end if                           
c
      if (limout) write(iwr,601)nelec,nconf,mstruc,mdetst,mstrco,mdetco,
     *                         lword,kend
601   format(/,5x,' parameters part 2 ...  ',/,
     *         5x,' nelec  ',i6,4x,' nconf  ',i6,4x,' mstruc ',i6,4x,
     *                                              ' mdetst ',i6,4x,/,
     *         5x,' mstrco ',i6,4x,' mdetco ',i6,4x,' lword  ',i10,4x,
     *            '  kend  ',i10)
c
      if (kend.gt.lword) call corfait(kend,lword,'crestr before crestd')
c
c... clean irumer array (called q(kbonds) here)
c
      call izero(nconf*maxsp*nelec,q(kbonds),1)
c
      if (mrsd) then
         call mrsdd (q(kconf),melec1,nsp0,q(nspath),q(ispath),
     &               maxsp,q(korb),ninner,q(kvirt),nvirt,q(idconf),
     &               q(knstru),q(kntdet),q(kdetps),q(idetps),q(kcofff),
     &               q(kpacde),q(kstruc),q(kcoeff),q(kdeter),q(klogic),
     &               q(kparit),q(kkscr),mscr,q(kscr),mconf,q(kbonds))
      else
         call crestd(q(kconf) ,q(knstru),q(kntdet),q(kstruc),
     *               q(kcoeff),q(kdeter),q(klogic),q(nspath),
     *               q(ispath),maxsp,q(kcofff),q(kdetps),q(idetps),
     *               q(kpacde),q(kkscr),mscr,q(kscr),melec1,nconf,ntdet,
     *               q(kbonds))
      end if
_IF(atmol)
      close(15)
_ENDIF
c
      return
      end
c***********************************************************************
      subroutine crterm (kterm,msocc,icount,msp)
c
      implicit REAL  (a-h,o-z) , integer (i-n)
c
INCLUDE(common/turtleparam)
c
c***********************************************************************
c
c     'crterm' generates all leading terms that can be generated with
c     'nalpha' alpha-electrons and 'nbeta' beta-electrons. (these
c     numbers are given through by the calling program by means of com-
c     mon-block 'condep').
c     the created leading terms are stored as columns of the array
c     'kterm', each row-element being '1' or '0', standing for alpha-
c     and beta-spin respectively.
c
c     for the scheme in which the leading terms are generated see
c           m.raimondi,m.simonetta,g.f.tantardini,
c                        'ab initio valence bond theory'
c      computer physics report 2 (1981), p192-194
c
c     msocc === nsocc        )
c
      dimension kterm (msocc, * )
c
      common /condep/ nsocc,ndocc,nalpha,nbeta,npair
c
c            ****   alpha-spin :  1   ****
c            ****    beta-spin :  0   ****
c
c                       create first leading term
c
      do 10  iodd=1,nsocc,2
         kterm(iodd,1) = 1
10    continue
c
      do 20 ieven=2,nsocc,2
         if (ieven.le.nbeta*2) then
            kterm(ieven,1) = 0
         else
            kterm(ieven,1) = 1
         endif
20    continue
c
      icount = 1
c
c                    create all other leading terms
c     (construction of a new term, is dependent of the previous term)
c
c      write(iwr,*)'remco nsocc',nsocc,nbeta,msp
30    if (icount.gt.maxspvb) call vberr('maxspvb(turtleparam) exceeded')
c
c                       find lowest beta-element
c
      jcheck = 2
c
40    continue
      if (kterm(jcheck,icount).ne.0) then
c
c                   this is no beta, try next element
c
         jcheck = jcheck + 1
         goto 40
      else
c
         nseqb = 1
c
c         jump out of routine if last leading term was already made
c.ab      jump out of routine if highest needed leading term
c.ab        was already created
c
         if ((jcheck.eq.nsocc-nbeta+1).or.(icount.eq.msp)) goto 100
c
c     find first alpha-element, that is preceded by a beta-element
c
50       continue
         if (kterm(jcheck+1,icount).ne.1) then
c
c       the sequence of beta's appears to contain another element
c
            nseqb = nseqb + 1
c
c                          try next element
c
            jcheck = jcheck + 1
            goto 50
         else
c
            if (nseqb.lt.2) then
c
c         copy elements 1 to jcheck-1 from previous leading term
c
               do 60 icp=1,jcheck-1
                  kterm(icp,icount+1) = kterm(icp,icount)
60             continue
            else
c
c     put beta's 1 to nseqb-1 on their start-positions (that is: as in
c     the first term), other elements before jcheck-1 must be alpha's
c
               do 70  iodd=1,jcheck-1,2
                  kterm(iodd,icount+1) = 1
70             continue
c
               do 80 ieven=2,jcheck-1,2
                  if (ieven.le.(nseqb-1)*2) then
                     kterm(ieven,icount+1) = 0
                  else
                     kterm(ieven,icount+1) = 1
                  endif
80             continue
            endif
c
c     copy the elements jcheck and jcheck+1 exchanged from previous term
c
            kterm(jcheck  ,icount+1) = kterm(jcheck+1,icount)
            kterm(jcheck+1,icount+1) = kterm(jcheck  ,icount)
c
c       copy the elements jcheck+2 to nsocc from previous leading term
c
            do 90 icp=jcheck+2,nsocc
               kterm(icp,icount+1) = kterm(icp,icount)
90          continue
c
         endif
      endif
c
c        a leading term has just been created, proceed with next one
c
      icount = icount + 1
      goto 30
c
100   continue
c
      return
      end
      subroutine deter (kposar,ndetst,kbeta,kconf,knwdet,nnwdet,
     *                  nelec,npppr,mperm,nnwstr,knwstr,mdetst,
     *                  klogic)
c***********************************************************************
c
      implicit REAL  (a-h,o-z) , integer (i-n)
c
c***********************************************************************
c
c     'deter' creates all determinants, that can be generated by assig-
c     nign beta-spin to 'nbeta' orbitals out of 'nsocc' single occupied
c     orbitals. for each determinant, with the help of the array
c     'kposar', a reference to the array 'kbeta' is be calculated. this
c     array contains information about which orbitals of the configura-
c     tion (which is defined by the array 'kconf)
c     must be give beta-spin.
c     npppr = npair
c
      dimension kposar ( * ),
     &          kbeta  ( npppr  , mperm , 0:npppr   ,  *  ),
     &          kconf  ( nelec+1 ),
     &          knwdet ( nelec   ,  *  ),
     &          knwstr ( mdetst,   *   ),
     &          klogic ( * )
c
INCLUDE(common/turtleparam)
      common /   fac/ kfacul(0:mmsocc)
      common /condep/ nsocc,ndocc,nalpha,nbeta,npair
c
ckoos      nnwdet = kfacul(nsocc)/(kfacul(nsocc-nbeta)*kfacul(nbeta))
c.ab      nnwdet = newtont(nsocc,nbeta)
      nnwdet = 0
ckoos
      nperm  = nbeta
c
      do 10 istruc=1,nnwstr
      do 15 idet=1,ndetst
c
         nnwdet = nnwdet + 1
         ndtemp = locat1(klogic,nnwdet-1,knwstr(idet,istruc))
         if (ndtemp.eq.0) then
            klogic(nnwdet) = knwstr(idet,istruc)
            jposar = kposar(klogic(nnwdet))
c
c     determine the last three indices of the array kbeta, wherewith the
c     beta-orbitals belonging to the current determinant can be found
c
c                             fourth index
c
         index4 = 1
c
20       continue
         if (.not.(jposar.gt.(index4-1)*ndetst
     &       .and.jposar.le.(index4  )*ndetst)) then
            index4 = index4 + 1
            goto 20
         endif
c
c                             third index
c
          jcomp = jposar - (index4-1)*ndetst
         index3 = 0
         jlower = 0
         jupper = 1
c
30       continue
         if (.not.(jcomp.gt.jlower
     &       .and.jcomp.le.jupper)) then
            index3 = index3 + 1
            jlower = jupper
ckoos            jupper = jupper + kfacul(nperm)/
ckoos     &                     (kfacul(nperm-index3)*kfacul(index3))
            jupper = jupper + newtont(nperm,index3)
ckoos
            goto 30
         endif
c
c                             second index
c
         index2 = jcomp - jlower
c
c                          define determinant
c
c                       single occupied orbitals
c
         jalpha = 1
          jbeta = 1
c
         do 40 isocc=1,nsocc
            if ((isocc.ne.kbeta(jbeta,index2,index3,index4)).or.
     *          (jbeta.gt.nbeta)) then
c
c                           this is an alpha
c
               knwdet(             jalpha,nnwdet) = kconf(isocc)
               jalpha = jalpha + 1
            else
c
c                           this is a beta
c
               knwdet(nalpha+ndocc+jbeta ,nnwdet) = kconf(isocc)
                jbeta =  jbeta + 1
            endif
40       continue
c
c                      double occupied orbitals
c
         do 50 idocc=1,ndocc
c
c     put respectively an alpha- and an (identical) beta-orbital on the
c                                                              array
c
            knwdet(nalpha            +idocc,nnwdet)=
     &                                       kconf(nsocc+2*idocc)
            knwdet(nalpha+ndocc+nbeta+idocc,nnwdet)=
     &                                       kconf(nsocc+2*idocc)
50       continue
c.ab
         else
           nnwdet = nnwdet - 1
         endif
c
15    continue
10    continue
c
      return
      end
      subroutine dimerr(string,number)
c...  quit because of dimension error
      implicit REAL (a-h,o-z), integer (i-n)
c
      character*(*) string
INCLUDE(../m4/common/iofile)
c
      write(iwr,5)
5     format(//30x,'error!!!!!!',//10x,'dimensions are overwritten')
      write(iwr,15) string,number 
15    format(//'the number of ',a,' exceeds',i6)
c
      stop
      end
      subroutine doubc(kconf ,melec1,nrefco,nspath,ispath,maxsp , 
     &                 korb  ,norb  ,idconf,
     &                 jact1 ,jvirt1,jact2 ,jvirt2,
     &                 knstru,kntdet,kdetps,idetps,coefff,pacdet,
     &                 nconf ,nstruc,ncfff ,npacd ,ndet  ,
     &                 kstruc,coeff ,kdeter,klogic,parity,kscr,
     &                 kcurcf,mscr,scr   ,mconf,irumer)
c...
c...                        jact1 --> jvirt1
c...                        jact2 --> jvirt2
c...  'doubc' performs the hereabove illustrated double excitation to
c...  reference configurations. the doubles are not stored, but spin 
c...  projected immediately after their creation. 
c...
c...  the reference configurations and the (temporarily stored) excited 
c...  configurations are constructed as follows:
c...       /        singles block        /        doubles block        /
c...               unique numbers                pairs of numbers
c...  size:           nsocc                         2*ndocc
c...
c...  this leads to the following construction of determinants:
c...       /        alpha block        /         beta block        /
c...       /   singles   /   doubles   /   singles   /   doubles   /
c...  size: nalpha-ndocc     ndocc      nbeta-ndocc      ndocc
c...  (both doubles blocks are the same)
c...
c...  within each block, orb(i+1) >= orb(i) holds (for both  
c...  configurations and determinants). this ordering can only be
c...  disturbed in determinants by excitation to external virtuals, for
c...  these are always put at the last positions of their spin blocks.
c...  
      implicit REAL (a-h,o-z), integer (i-n)
c
INCLUDE(common/c8_16vb)
INCLUDE(common/turtleparam)
c
      common /indep/ time0,nelec,mult
c                                        
INCLUDE(common/logics)
c
      dimension kconf (melec1,*)
      dimension nspath(*),
     &          ispath(maxsp,*),
c
     &          korb  (*),
     &          idconf(*),
c
     &          kactiv(2),
     &          kvirt (2),
c
     &          knstru(*),
     &          kntdet(*),
     &          kdetps(*),
     &          idetps(*),
     &          coefff(*),
     &          pacdet(*),
c
     &          kstruc(*),
     &          coeff (*),
     &          parity(*),
     &          kdeter(melec1-1,*),
     &          klogic(*),
     &          kcurcf(*),
     &          kscr  (*),
     &          irumer(*),
     &          scr   (*)
c
      logical nofail
c
      nalpha = (nelec + mult - 1)/2
      nofail = .false.
c
      kactiv(1) = jact1
      kactiv(2) = jact2
      kvirt(1)  = jvirt1
      kvirt(2)  = jvirt2
c
      do 10 irefco=1,nrefco
c
         ndocc  = kconf(nelec+1,irefco)
         nsocc  = nelec - 2*ndocc
c...  determine current positions of actives and virtuals
         jposa1 = ifind(kconf(1,irefco),nelec,kactiv(1))
         jposa2 = ifind(kconf(1,irefco),nelec,kactiv(2))
         jposv1 = ifind(kconf(1,irefco),nelec,kvirt(1))
         jposv2 = ifind(kconf(1,irefco),nelec,kvirt(2))
c...  leave out excitations that lead to nul-functions
         if ( jposa1.eq.0.or.
     &        jposa2.eq.0.or.
     &       (jposa1.eq.jposa2.and.jposa1.le.nsocc).or.
     &        jposv1.gt.nsocc.or.
     &        jposv2.gt.nsocc.or.
     &       (jposv1.eq.jposv2.and.jposv1.gt.0))         go to 10
         call icopy(nelec+1,kconf(1,irefco),1,kcurcf,1)
c...  perform two excitations sequentially        
         do 15 iexc=1,2
            jactiv = kactiv(iexc)
            jvirt  = kvirt (iexc)
            if (iexc.eq.1) then
               jposa = jposa1
               jposv = jposv1
            else
               jposa = ifind (kcurcf,nelec,jactiv)
               jposv = ifindo(kcurcf,nsocc,jvirt)  
            end if
            call excitc(jactiv,jvirt,jposa,jposv,nsocc,ndocc,kcurcf)
15       continue
c...  define structures from new configuration
         call projc(kcurcf,melec1,1,irefco,nspath(irefco),
     &              ispath(1,irefco),maxsp,nofail,knstru(nconf+1),
     &              kntdet(nconf+1),kdetps(nstruc+1),idetps(ncfff+1),
     &              coefff(ncfff+1),pacdet(npacd+1),nnstr,nncfff,nnpac,
     &              nndet,krejcf,nrejcf,kstruc,coeff,kdeter,klogic,
     &              kscr,mscr,scr,irumer)
         if (nrejcf.eq.1) go to 10
c
         if ((ifind(kvirt,2,2**n8_16-2).gt.0)
     &       .and.(kvirt(1).ne.kvirt(2))) then
c...  there are one or more external virtuals in the singles parts of an  
c...  alpha and/or beta block. put the externals last in their spin  
c...  blocks. (disturb the orbital-ordering within the determinants). 
c...  in case of jact1 --> 254, jact2 --> 254, externals have been put
c...  last in their spin blocks automatically
            call unpack(pacdet(npacd+1),n8_16,kdeter,
     &                  kntdet(nconf+1)*nelec)
            if (kvirt(2).eq.2**n8_16-1) then
c...  jact1 --> 254
c...  jact2 --> 255
               next = 2
            else
c...  jact1 --> internal
c...  jact2 --> 254
               next = 1
            end if
            do 110 idet=1,kntdet(nconf+1)
               nperm = 0
               do 120 iexc=1,next
                  if (iexc.eq.1) then
                     jvirt = 2**n8_16-2
                  else
                     jvirt = 2**n8_16-1
                  end if
                  jpose = ifind(kdeter(1,idet),nelec,jvirt)
                  if (jpose.le.nalpha) then
                     npose = nalpha
                  else
                     npose = nelec 
                  end if
                  call open1(kdeter(1,idet),npose,jpose)
                  kdeter(npose,idet) = jvirt
                  nperm = nperm + iabs(npose-jpose)
120            continue
               parity(idet) = (-1.0d0)**nperm
110         continue
            call pack(pacdet(npacd+1),n8_16,kdeter,
     &                kntdet(nconf+1)*nelec)
            do 130 icfff=ncfff+1,ncfff+nncfff
               coefff(icfff) = parity(idetps(icfff))*coefff(icfff)
130         continue
         end if
c
c...  raise configuration counters 
         nconf  = nconf  + 1
         if (nconf.eq.mconf) call dimerr('configurations',mconf)
         ncfff  = ncfff  + nncfff
         nstruc = nstruc + nnstr
         npacd  = npacd  + nnpac
         ndet   = ndet   + nndet
c...  store configuration identification
         call storid(kcurcf,nconf,idconf,korb,norb)
c
10    continue
c
      return
      end
      subroutine doubex (kstate,ndim,norb,jcolgr,jcolnw,kexcit,nexcit,
     &                                                         ngener)
c***********************************************************************
c
      implicit REAL  (a-h,o-z) , integer (i-n)
c
c***********************************************************************
c
      dimension kstate (ndim, * ),
     &          kexcit ( * )
c
      jstate =                jcolnw
      nexref =      kexcit(nexcit+1)
       ndocc = kstate(norb+1,jcolgr)
       nsocc =        norb - 2*ndocc
c
      do 10 iskip1=1,norb
c*****
c
c     the double occupied orbitals are stored twice, after one another
c
         if (iskip1.gt.nsocc) then
            if (((iskip1-nsocc)/2)*2.eq.iskip1-nsocc) goto 5
         endif
c
         do 20 iskip2=iskip1+1,norb
c********
c
c     the double occupied orbitals are stored twice, after one another
c
            if (iskip2.gt.nsocc) then
               if (((iskip2-nsocc)/2)*2.ne.iskip2-nsocc) goto 15
            endif
c
            do 30 iexc1=1,nexcit
c***********
               do 40 iexc2=iexc1,nexcit
c**************
                  if (iexc1.gt.nexcit-nexref) then
c
c     no more than     no more than one excitation is possible to a virt
c
                     if (iexc1.eq.iexc2) goto 25
c
c                  an orbital cannot be excited to itself
c
                     if (iskip1-(nsocc-nexref).eq.iexc1-(nexcit-nexref))
     &                                                           goto 25
                     if (iskip1-(nsocc-nexref).eq.iexc2-(nexcit-nexref))
     &                                                           goto 25
                     if (iskip2-(nsocc-nexref).eq.iexc1-(nexcit-nexref))
     &                                                           goto 25
                     if (iskip2-(nsocc-nexref).eq.iexc2-(nexcit-nexref))
     &                                                           goto 25
                  endif
c
                  if (iexc2.gt.nexcit-nexref) then
c
c                  an orbital cannot be excited to itself
c
                     if (iskip1-(nsocc-nexref).eq.iexc2-(nexcit-nexref))
     &                                                           goto 25
                     if (iskip2-(nsocc-nexref).eq.iexc2-(nexcit-nexref))
     &                                                           goto 25
                  endif
c
c                      perform double excitations
c
                  do 50 icp=1,iskip1-1
                     kstate(icp,jstate) = kstate(icp,jcolgr)
50                continue
                  kstate(iskip1,jstate) = kexcit(iexc1)
                  do 60 icp=iskip1+1,iskip2-1
                     kstate(icp,jstate) = kstate(icp,jcolgr)
60                continue
                  kstate(iskip2,jstate) = kexcit(iexc2)
                  do 70 icp=iskip2+1,norb
                     kstate(icp,jstate) = kstate(icp,jcolgr)
70                continue
c
c     reorder configuration by putting the double occupied orbitals af-
c     ter the single occupieds. store the number of double occupied or-
c     bitals as row-element 'norb + 1'
c
                  call occup (kstate(1,jstate),norb)
c
                  jstate = jstate + 1
c
25                continue
40             continue
30          continue
15          continue
20       continue
5        continue
10    continue
c
      ngener = jstate - jcolnw
c
      return
      end
      subroutine doubs(kconf ,melec1,nrefco,
     &                 korb  ,norb  ,idconf,
     &                 jact1 ,jvirt1,jact2 ,jvirt2,
     &                 knstru,kntdet,kdetps,idetps,coefff,pacdet,
     &                 nconf ,nstruc,ncfff ,npacd ,ndet  ,
     &                 kdeter,parity,krejd ,keqdet,kscr  ,scr   ,kcurcf,
     &                 mscr  ,mconf)
c...                                                                   
c...  'doubs' performs the double excitation jact1 --> jvirt1,
c...  jact2 --> jvirt2 to reference structures, by applying the 
c...  following excitation operator twice:   
c...            e(a-->v) = v(a)+a(a) + v(b)+a(b)
c...  (v=virtual, a=active, (a)=alpha, (b)=beta, +=cross)
c...
c...  assumed reference configuration-construction:
c...       /        singles block        /        doubles block        /
c...               unique numbers                pairs of numbers
c...  size:           nsocc                         2*ndocc  
c...  within each block, orb(i+1) > orb(i) holds for all i
c...
c...  assumed reference determinant-construction:
c...       /        alpha block        /         beta block        /
c...       /   singles   /   doubles   /   singles   /   doubles   /
c...  size:  nalpha-ndocc     ndocc      nbeta-ndocc      ndocc
c...  (both doubles blocks are the same)
c...  within each block, orb(i+1) > orb(i) holds for all i. this
c...  ordering might be disturbed if external virtuals are involved. 
c...
c...  external virtuals are put at the last positions in their spin 
c...  blocks
c...
      implicit REAL (a-h,o-z), integer (i-n)
c 
INCLUDE(common/c8_16vb)    
c
      common /indep/ time0,nelec,mult
c
      logical null
c
      dimension kconf (melec1,*),
c
     &          korb  (*),
     &          idconf(*),
     &          kactiv(2),
     &          kvirt (2),

c
     &          knstru(*),
     &          kntdet(*),
     &          kdetps(*),
     &          idetps(*),
     &          coefff(*),
     &          pacdet(*),
c
     &          kdeter(melec1-1,*),
     &          parity(*),
     &          krejd (*),
     &          keqdet(*),
     &          kscr  (*),    
     &          scr   (*),
     &          kcurcf(*)
c
      kactiv(1) = jact1
      kactiv(2) = jact2
      kvirt(1)  = jvirt1
      kvirt(2)  = jvirt2
      nalpha = (nelec + mult - 1)/2
      nbeta  = nelec - nalpha
c...  initialize adresses base reference 
      jstruc = 1
      jcfff  = 1
      jpacd  = 1
c...  for the second excitation single excited structures are used as
c...  reference!!!
      do 10 irefco=1,nrefco
         ndocc = kconf(nelec+1,irefco)
         nsocc = nelec - 2*ndocc
c...  determine current positions actives and virtuals
         jposa1 = ifind(kconf(1,irefco),nelec,kactiv(1))
         jposa2 = ifind(kconf(1,irefco),nelec,kactiv(2))
         jposv1 = ifind(kconf(1,irefco),nelec,kvirt(1) )
         jposv2 = ifind(kconf(1,irefco),nelec,kvirt(2) )
c...  leave out excitations leading to nul-functions
         if (jposa1.eq.0.or.jposa2.eq.0.or.
     &       (jposa1.eq.jposa2.and.jposa1.le.nsocc).or.
     &       jposv1.gt.nsocc.or.jposv2.gt.nsocc.or.
     &       (jposv1.eq.jposv2.and.jposv1.gt.0)) go to 10000
         call icopy(nelec+1,kconf(1,irefco),1,kcurcf,1)
         call unpack(pacdet(jpacd),n8_16,kdeter,nelec*kntdet(irefco))
c...  loop over excitations
         do 15 iexc=1,2
            jactiv = kactiv(iexc)
            jvirt  = kvirt (iexc)
            if (iexc.eq.1) then
               jposa = jposa1
               jposv = jposv1
c...  reference adresses
               jrconf = irefco
               jrstru = jstruc
               jrcfff = jcfff
c...  size of reference
               nrcfff = 0
               do 20 istruc=jrstru,jrstru-1+knstru(irefco)
                  nrcfff = nrcfff + kdetps(istruc)
20             continue
c...  temporal storage of singly excited structures
               jnconf = nconf + 2
               jnstru = nstruc +   knstru(irefco) + 1
               jncfff = ncfff  + 4*nrcfff + 1
            else
               jposa = ifind (kcurcf,nelec,jactiv)
               jposv = ifindo(kcurcf,nsocc,jvirt)  
c...  reference adresses (=singly excited structures)
               jrconf = jnconf
               jrstru = jnstru
               jrcfff = jncfff
c...  size of reference
               nrcfff = nncfff
c...  start adresses doubly excited structures
               jnconf = nconf  + 1
               jnstru = nstruc + 1
               jncfff = ncfff  + 1
            end if
c...  go for it
            call excits(jactiv,jvirt,jposa,jposv,nalpha,nbeta,ndocc,
     &                  nsocc,jrconf,jrstru,jrcfff,nrcfff,jnconf,jnstru,
     &                  jncfff,nncfff,knstru,kntdet,kdetps,idetps,
     &                  coefff,null,kdeter,melec1,parity,krejd,keqdet,
     &                  kscr,scr,mscr)
            if (null) go to 10000
            call reconf(kcurcf,melec1,nelec,kdeter(1,1),1)
            call ordernumbs(kcurcf,nsocc)
15       continue
         if ((ifind(kvirt,2,2**n8_16-2).gt.0.or.
     &        ifind(kvirt,2,2**n8_16-1).gt.0)
     &        .and.(kvirt(1).ne.kvirt(2))        ) then
c...  there are one or more external virtuals in the singles part of an  
c...  alpha or beta block. now the ordering of the orbitals is disturbed 
c...  to put the externals last in their spin blocks. 
            if (ifind(kvirt,2,2**n8_16-1).eq.0) then
               next  = 1
               jvirt = 2**n8_16-2
            else
               next  = 2
            end if
            do 1010 idet=1,kntdet(jnconf)
               nperm = 0
               do 1020 iexc=1,next
                  if (next.eq.2) then
                     if (iexc.eq.1) then
                        jvirt = 2**n8_16-2
                     else
                        jvirt = 2**n8_16-1
                     end if
                  end if
                  jpose = ifind(kdeter(1,idet),nelec,jvirt)
                  if (jpose.le.nalpha) then
                     npose = nalpha
                  else
                     npose = nelec
                  end if
                  call open1(kdeter(1,idet),npose,jpose)
                  kdeter(npose,idet) = jvirt
                  nperm = nperm + iabs(npose-jpose)
1020           continue
               parity(idet) = (-1.0d0)**nperm
1010        continue
            do 1030 icfff=ncfff+1,ncfff+nncfff
               coefff(icfff) = parity(idetps(icfff))*coefff(icfff)
1030        continue
         end if
c...  write determinants
         call pack(pacdet(npacd+1),n8_16,kdeter,nelec*kntdet(jnconf))
c...  normalize structures
         call normstruc(kdetps(jnstru),coefff(jncfff),knstru(jnconf))
c...  update size parameters
         nconf  = nconf  + 1
         if (nconf.eq.mconf) call dimerr('configurations',mconf)
         nstruc = nstruc + knstru(jnconf)
         ncfff  = ncfff  + nncfff
         npacd  = npacd  + (nelec*kntdet(jnconf)-1)/(64/n8_16) + 1
         ndet   = ndet   + kntdet(jnconf)
c...  store configuration specifications
         call storid(kcurcf,nconf,idconf,korb,norb)
10000    continue
c....  calculate once again size of base reference 
         nrcfff = 0
         do 10020 istruc=jstruc,jstruc-1+knstru(irefco)
            nrcfff = nrcfff + kdetps(istruc)
10020    continue
c...  raise adresses base reference 
         jstruc = jstruc + knstru(irefco)
         jcfff  = jcfff  + nrcfff
         jpacd  = jpacd  + (nelec*kntdet(irefco)-1)/(64/n8_16) + 1
10    continue
c
      return
      end    
      subroutine eqdet(ksbeta,nsbeta,ndet,krejd,nprejd,keqdet,maxeq,
     &                 nequal)
c...  completely devoted to excits
c...    
c...  'eqdet' gathers groups of equal determinants on keqdet. 
c...  the comparance of the determinants is restricted to the single
c...  occupied beta orbitals, gathered on ksbeta. it is supposed that
c...  bsbeta(i+1,idet) > ksbeta(i,idet) for all i.
c...
c...  in the first part keqdet is given the following structure:
c...  the number of elements in a group is followed by its elements,
c...  which are followed again by the number of elements of a new group 
c...  and so on, and so on. 
c...  
c...  in the second part keqdet is rearranged to get an essentially 
c...  2-dimensional structure:
c...                               keqdet(1,ieq): number of elements
c...                                                        in group
c...  keqdet(2,ieq)... keqdet(keqdet(1,ieq),ieq): equal determinants
c...                       keqdet(1+maxeq+1,ieq): 0 (has a purpose in
c...                                                   singst/doubst)
c...
c...  the pauli-forbidden determinants are not added to ksbeta, but
c...  form however members of kdeter. to account for this shift-values
c...  are introduced to relate the two types of adresses:  
c...  idet1/2:            adress for ksbeta 
c...  idet1/2 + jshif1/2: adress for kdeter (in singst/doubst)
c...  (ndet is already updated)
c...
      implicit REAL (a-h,o-z), integer (i-n)
c
      dimension ksbeta(nsbeta,*),
     &          keqdet(*),
     &          krejd (*)
c...
c...                      group determinants
c...
      nequal = 0
      maxeq  = 0
      nextgr = 1
c
      jprej1 = 1
      jshif1 = 0
      do 10 idet1=1,ndet-1
c...  calculate adress shift for idet1 
20       if (jprej1.le.nprejd) then
            if (idet1+jshif1.eq.krejd(jprej1)) then
               jshif1 = jshif1 + 1
               jprej1 = jprej1 + 1 
               go to 20
            end if
         end if
c...  check if idet1 isn't already in a group
         jjj = 1
         do 30 iequal=1,nequal
            do 40 ieqdet=1,keqdet(jjj)
               if (idet1+jshif1.eq.keqdet(jjj+ieqdet)) go to 10
40          continue
            jjj = jjj + keqdet(jjj) + 1  
30       continue
c
         keqdet(nextgr) = 0
         jprej2 = jprej1
         jshif2 = jshif1
         do 50 idet2=idet1+1,ndet
c...  calculate adress shift for idet2 
60          if (jprej2.le.nprejd) then
               if (idet2+jshif2.eq.krejd(jprej2)) then
                  jshif2 = jshif2 + 1
                  jprej2 = jprej2 + 1 
                  go to 60
               end if
            end if
c...  check if idet2 isn't already in a group
            jjj = 1
            do 70 iequal=1,nequal
               do 80 ieqdet=1,keqdet(jjj)
                  if (idet2+jshif2.eq.keqdet(jjj+ieqdet)) go to 50
80             continue
               jjj = jjj + keqdet(jjj) + 1 
70          continue
c...  compare determinants
            do 90 iorb=1,nsbeta
               if (ksbeta(iorb,idet1).ne.ksbeta(iorb,idet2)) go to 50
90          continue
c...  idet1 equals idet2 
            if (keqdet(nextgr).eq.0) then 
c...  this is the beginning of a new group
               keqdet(nextgr+1) = idet1 + jshif1
               keqdet(nextgr  ) = 1 
               nequal = nequal + 1
            end if
            keqdet(nextgr+keqdet(nextgr)+1) = idet2 + jshif2
            keqdet(nextgr) = keqdet(nextgr) + 1
c
50       continue
         if (keqdet(nextgr).gt.0) then
            maxeq  = max(maxeq,keqdet(nextgr))
            nextgr = nextgr + keqdet(nextgr) + 1 
         end if
10    continue
c...
c...      turn keqdet into an effectively 2-dimensional array
c...
      do 110 iequal=nequal,1,-1
c...  calculate current position of group iequal  
         jcurr = 1
         do 120 i=1,iequal-1
            jcurr = jcurr + keqdet(jcurr) + 1 
120      continue
c...  calculate new position of group iequal
         jnew = (iequal-1)*(maxeq+2) + 1
c...  move group iequal
         call copi(keqdet(jcurr),keqdet(jnew),1+keqdet(jcurr),-1)
         keqdet(jnew+keqdet(jnew)+1) = 0
110   continue
c
      return
      end
      subroutine excitc(jactiv,jvirt,jposa,jposv,nsocc,ndocc,kcurcf)
c...  
c...  jactiv: active
c...   jvirt: virtual
c...   jposa: position of active  on kcurcf
c...   jposv: position of virtual on kcurcf
c...   nsocc: number of single occupieds
c...   ndocc: number of double occupieds
c...  kcurcf: configuration
c...
c...  'excitc' performs the excitation jactiv --> jvirt on kcurcf.
c...  the structure of kcurcf is as follows:
c...       /        singles block        /        doubles block        /
c...               unique numbers                pairs of numbers
c...  size:           nsocc                         2*ndocc
c...
c...  within each block, orb(i+1) >= orb(i) holds,and this remains so 
c...  after the excitation. 
c...  nsocc and ndocc are updated.
c...  
      implicit REAL (a-h,o-z), integer (i-n)
c
      common /indep/ time0,nelec,mult
c
      dimension kcurcf(*)
c...
c...  select proper excitation type according to degree of occupation 
c...  active and virtual
c...
      if (jposv.eq.0) then
         if (jposa.le.nsocc) then 
            jcase = 1
         else
            jcase = 2
         end if
      else
         if (jposa.le.nsocc) then 
            jcase = 3
         else
            jcase = 4
         end if
      end if
      go to (1,2,3,4) jcase
c
1     continue
c...
c...  jcase 1: virtual is not occupied, active is singly occupied
c...
      nposv = irank(kcurcf,nsocc,jvirt,jposa,0)
      call open1(kcurcf,nposv,jposa)
      kcurcf(nposv) = jvirt
      return
c
2     continue
c...
c...  jcase 2: virtual is not occupied, active is doubly occupied
c...
      nposa = irank(kcurcf,nsocc,jactiv,0,0)
      nposv = irank(kcurcf,nsocc,jvirt ,0,0)
      if (jactiv.lt.jvirt) then
         nposv = nposv + 1
         call open1(kcurcf,nposa,jposa)
         kcurcf(nposa) = jactiv
         call open1(kcurcf,nposv,jposa+1)
         kcurcf(nposv) = jvirt
      else
         nposa = nposa + 1
         call open1(kcurcf,nposv,jposa)
         kcurcf(nposv) = jvirt
         call open1(kcurcf,nposa,jposa+1)
         kcurcf(nposa) = jactiv
      end if
      ndocc = ndocc - 1
      nsocc = nelec - 2*ndocc
      kcurcf(nelec+1) = ndocc
      return
c
3     continue
c...
c...  jcase 3: virtual is singly occupied, active is singly occupied
c...
      nposv = nsocc - 2 + irank(kcurcf(nsocc+1),2*ndocc,jvirt,0,0)
      call open1(kcurcf,nposv+1,max(jposa,jposv))
      kcurcf(nposv+1) = jvirt
      call open1(kcurcf,nposv  ,min(jposa,jposv))
      kcurcf(nposv)   = jvirt
      ndocc = ndocc + 1
      nsocc = nelec - 2*ndocc
      kcurcf(nelec+1) = ndocc
      return
c
4     continue
c...
c...  jcase 4: virtual is singly occupied, active is doubly occupied
c...
      nposa = irank(kcurcf,nsocc,jactiv,jposv,0)
      call open1(kcurcf,nposa,jposv)
      kcurcf(nposa) = jactiv
      nposv = nsocc + irank(kcurcf(nsocc+1),2*ndocc,jvirt,
     &                                  jposa-nsocc,jposa-nsocc+1)
      if (nposv.le.jposa) then
         call open1(kcurcf,nposv  ,jposa)
         kcurcf(nposv)   = jvirt
         call open1(kcurcf,nposv+1,jposa+1)
         kcurcf(nposv+1) = jvirt
      else
         call open1(kcurcf,nposv+1,jposa+1)
         kcurcf(nposv+1) = jvirt
         call open1(kcurcf,nposv  ,jposa)
         kcurcf(nposv)   = jvirt
      end if
      return
c
      end
      subroutine excits(jactiv,jvirt ,jposa ,jposv ,
     &                  nalpha,nbeta ,ndocc ,nsocc ,
     &                  jrconf,jrstru,jrcfff,nrcfff,
     &                  jnconf,jnstru,jncfff,nncfff,
     &                  knstru,kntdet,kdetps,idetps,coefff,
     &                  null  ,
     &                  kdeter,melec1,parity,krejd ,keqdet,kscr  ,scr  ,
     &                  mscr)
c...                                                                   
c...  'excits' applies the following excitation operator to reference 
c...  structures:  e(a-->v) = v(a)+a(a) + v(b)+a(b)
c...  (v=virtual, a=active, (a)=alpha, (b)=beta, +=cross)
c...
c...  assumed reference determinant-construction:
c...       /        alpha block        /         beta block        /
c...       /   singles   /   doubles   /   singles   /   doubles   /
c...  size:  nalpha-ndocc     ndocc      nbeta-ndocc      ndocc
c...  (both doubles blocks are the same)
c...  within each block, orb(i+1) > orb(i) holds. this ordering is 
c...  maintained after the excitation; the number of permutations needed
c...  for this maintainance might cause a parity change. 
c...
c...  start adresses of reference structures: jrconf, jrstru, jrcfff 
c...  size of reference structures: knstru(jrconf),kntdet(jrconf),nrcfff
c...  start adresses of generated structures: jnconf, jnstru, jncfff
c...  size of reference structures: knstru(jnconf),kntdet(jnconf),nncfff
c...
      implicit REAL (a-h,o-z), integer (i-n)    
c
      parameter (nav = 2)
c
      common /indep/ time0,nelec,mult
c
      logical null
c
      dimension knstru(*),
     &          kntdet(*),
     &          kdetps(*),
     &          idetps(*),
     &          coefff(*),
c
     &          kdeter(melec1-1,*),
     &          parity(*),
     &          krejd (*),
     &          keqdet(*),
     &          kscr  (*),
     &          scr   (*)
c
      null = .false.
c...  select proper excitation type according to degree of occupation 
c...  active and virtual
      if (jposv.eq.0) then
         if (jposa.le.nsocc) then 
            jcase = 1
         else
            jcase = 2
         end if
      else
         if (jposa.le.nsocc) then 
            jcase = 3
         else
            jcase = 4
         end if
      end if        
      go to (1,2,3,4) jcase
c...
c...        --------------------- excitation types ---------------------
c...
c
1        continue
c...
c...  jcase 1: virtual is not occupied, active is singly occupied
c...
c
c...
c...  copy and modify determinants
c...
c
c...  use educated guess for jposa at first determinant, based on
c...  ordering determinants after spin projection
      if (jposa.le.nbeta-ndocc) then
         jposa = nalpha + jposa
      else
         jposa = jposa - (nbeta-ndocc)
      end if
      do 100 idet=1,kntdet(jrconf)
c... search for position (jposa) of jactiv on column idet of kdeter
c... guess for jposa if idet exceeds 1: jposa(idet) = jposa(idet-1)
         call ifindd(kdeter(1,idet),jactiv,jposa,nalpha,ndocc)
c...  excite determinant
         if (jposa.le.nalpha) then
            nposv =          irank(kdeter(1       ,idet),
     &                             nalpha-ndocc,jvirt,
     &                             jposa      ,0)
         else
            nposv = nalpha + irank(kdeter(nalpha+1,idet),
     &                             nbeta -ndocc,jvirt,
     &                             jposa-nalpha,0)
         end if
         call open1(kdeter(1,idet),nposv,jposa)
         kdeter(nposv,idet) = jvirt
         parity(idet) = (-1.0d0)**(iabs(nposv-jposa))
100   continue
      kntdet(jnconf) = kntdet(jrconf)
c...
c...  copy structure-definitions
c...
      call icopy(nrcfff,idetps(jrcfff),1,idetps(jncfff),1)
      do 110 icfff=1,nrcfff
         coefff(jncfff-1+icfff) = 
     &   parity(idetps(jncfff-1+icfff))*coefff(jrcfff-1+icfff)
110   continue
      nncfff = nrcfff
      call icopy(knstru(jrconf),kdetps(jrstru),1,kdetps(jnstru),1)
      knstru(jnconf) = knstru(jrconf)
c
      return
c
2           continue
c...
c...  jcase 2: virtual is not occupied, active is doubly occupied
c...
c
c...
c...  copy and modify determinants
c...
      do 200 idet=kntdet(jrconf),1,-1
         do 210 iorb=1,nelec
            kdeter(iorb,2*idet-1) = kdeter(iorb,idet)
            kdeter(iorb,2*idet  ) = kdeter(iorb,idet)
210      continue
200   continue
c... jactiv is in doubles-part of both alpha- and beta-block, at
c... the same relative position. (the doubles-part is the same in all
c... determinants)
      jposa = nalpha - ndocc + (jposa-nsocc+1)/2
c...  replace alpha
      do 220 idet=1,kntdet(jrconf)
         nposv =          irank(kdeter(1       ,2*idet-1),
     &                          nalpha-ndocc,jvirt ,0,0)
         call open1(kdeter(1,2*idet-1),nposv,jposa)
         kdeter(nposv,2*idet-1) = jvirt
         nperm = iabs(nposv-jposa)
         nposa = nalpha + irank(kdeter(nalpha+1,2*idet-1),
     &                          nbeta -ndocc,jactiv,0,0)
         call open1(kdeter(1,2*idet-1),nposa,jposa+nbeta)
         kdeter(nposa,2*idet-1) = jactiv
         nperm = nperm + iabs(nposa-(jposa+nbeta))
         parity(2*idet-1) = (-1.0d0)**nperm
220   continue
c...  replace beta
      do 230 idet=1,kntdet(jrconf)
         nposv = nalpha + irank(kdeter(nalpha+1,2*idet)  ,
     &                          nbeta -ndocc,jvirt ,0,0)
         call open1(kdeter(1,2*idet)  ,nposv,jposa+nbeta)
         kdeter(nposv,2*idet)   = jvirt
         nperm = iabs(nposv-(jposa+nbeta))
         nposa =          irank(kdeter(1       ,2*idet)  ,
     &                          nalpha-ndocc,jactiv,0,0)
         call open1(kdeter(1,2*idet)  ,nposa,jposa)
         kdeter(nposa,2*idet)   = jactiv
         nperm = nperm + iabs(nposa-jposa)
         parity(2*idet) = (-1.0d0)**nperm
230   continue
      ndocc = ndocc - 1
      nsocc = nsocc + 2
      kntdet(jnconf) = 2*kntdet(jrconf)
c...
c...  copy and modify structure-definitions
c...
      do 240 icfff=1,nrcfff
         idetps(jncfff-1+2*icfff-1) = 2*idetps(jrcfff-1+icfff) - 1
         idetps(jncfff-1+2*icfff  ) = 2*idetps(jrcfff-1+icfff)
240   continue
      do 250 icfff=1,nrcfff
         coefff(jncfff-1+2*icfff-1) = 
     &   parity(idetps(jncfff-1+2*icfff-1))*coefff(jrcfff-1+icfff)
         coefff(jncfff-1+2*icfff  ) = 
     &   parity(idetps(jncfff-1+2*icfff  ))*coefff(jrcfff-1+icfff)
250   continue
      nncfff = 2*nrcfff
      do 260 istruc=1,knstru(jrconf)
         kdetps(jnstru-1+istruc) = 2*kdetps(jrstru-1+istruc)
260   continue
      knstru(jnconf) = knstru(jrconf)
c
      return
c
3           continue
c...
c...  jcase 3: virtual is singly occupied, active is singly occupied +++++
c...
c
c...
c...  copy and modify determinants
c...
      nprejd = 0
c...  use educated guess for jposa and jposv at first determinant, based
c...  on ordering determinants after spin projection
      if (jposa.le.nbeta-ndocc) then
         jposa = nalpha + jposa
      else
         jposa = jposa - (nbeta-ndocc)
      end if
      if (jposv.le.nbeta-ndocc) then
         jposv = nalpha + jposv
      else
         jposv = jposv - (nbeta-ndocc)
      end if
      do 300 idet=1,kntdet(jrconf)
c...  search for positions of jactiv (jposa) and jvirt (jposv) on column
c...  idet of kdeter
c...  guess for jposa/v if idet exceeds 1: jposa/v(idet) = 
c...                                       jposa/v(idet-1)
         call ifindd(kdeter(1,idet),jactiv,jposa,nalpha,ndocc)
         call ifindd(kdeter(1,idet),jvirt ,jposv,nalpha,ndocc)
c...  check wether jactiv and jvirt are in the same spin block
         if ((jposa.le.nalpha.and.jposv.le.nalpha).or.
     &       (jposa.gt.nalpha.and.jposv.gt.nalpha)) then
c...  the excitation is pauli-forbidden  
            krejd(nprejd+1) = idet
            nprejd = nprejd + 1
         else
c...  excite determinant
            nposv = nalpha - ndocc - 1 + irank
     &              (kdeter(nalpha-ndocc+1,idet),ndocc,jvirt,0,0)
            call open1(kdeter(1,idet),nposv      ,
     &                 min(jposa,jposv))
            kdeter(nposv      ,idet) = jvirt
            nperm = iabs(nposv-min(jposa,jposv))
            call open1(kdeter(1,idet),nposv+nbeta,
     &                 max(jposa,jposv))
            kdeter(nposv+nbeta,idet) = jvirt
            nperm = nperm + iabs((nposv+nbeta)-max(jposa,jposv))
            parity(idet) = (-1.0d0)**nperm
         end if
300   continue
      if (nprejd.eq.kntdet(jrconf)) then
         null = .true.
         return
      end if
      ndocc = ndocc + 1
      nsocc = nsocc - 2
c...  spin coupling of active and virtual in the reference structure
c...  have led to the multifold definition of determinants after the
c...  excitations. cope with this 
c...  gather single occupied beta orbitals
      nsbeta = nbeta - ndocc
      ii = 1
      iii = 1
      do 305 idet=1,kntdet(jrconf)
         if (iii.le.nprejd) then
            if (idet.eq.krejd(iii)) then
               iii = iii + 1
c...  forget this pauli-forbidden determinant
               go to 305
            end if
         end if
         call icopy(nsbeta,kdeter(nalpha+1,idet),1,kscr(ii),1)
         ii = ii + nsbeta
305   continue
c...  find equal determinants and group them
      call eqdet(kscr,nsbeta,kntdet(jrconf)-nprejd,krejd,nprejd,
     &           keqdet,maxeq,nequal)
c...  extend list of rejected (keep first of each equals group)
      ntrejd = nprejd
      do 310 iequal=1,nequal
         call icopy(keqdet((iequal-1)*(maxeq+2)+1)-1,
     &              keqdet((iequal-1)*(maxeq+2)+3),1,
     &              krejd(ntrejd+1),1)
         ntrejd = ntrejd + keqdet((iequal-1)*(maxeq+2)+1)-1
310   continue
      call ordernumbs(krejd(nprejd+1),ntrejd-nprejd)
      call icopy(ntrejd,krejd,1,kscr,1)
      call ordernumbs(kscr,ntrejd)
c...  remove rejected determinants
      if (ntrejd.gt.0)
     &   call skipi(kdeter,nelec,kntdet(jrconf),kscr,ntrejd,
     &              kdeter)
      kntdet(jnconf) = kntdet(jrconf) - ntrejd
c...
c...  copy and modify structure-definitions (mind rejected determinants)
c..
c
c...  renumber relevant determinants. output: kscr(oldnr) = newnr
      newnr = 1
      do 320 idet=1,kntdet(jrconf)
         if (ifindo(kscr,ntrejd,idet).eq.0) then
            kscr(ntrejd+idet) = newnr
            newnr = newnr + 1
         end if
320   continue
      call icopy(kntdet(jrconf),kscr(ntrejd+1),1,kscr(1),1)
c...  connect higher members of equals group with the first
      do 325 iequal=1,nequal
         jieq = (iequal-1)*(1+maxeq+1)
         nieq = keqdet(jieq+1)
         do 330  iieq=2,nieq
            kscr(keqdet(jieq+1+iieq)) = kscr(keqdet(jieq+1+1))
330      continue
325   continue
c...  store knstru, kdetps, idetps and coefff
      knstru(jnconf) = 0
      nncfff = 0
      jjcfff = jrcfff
      do 340 istruc=jrstru,jrstru-1+knstru(jrconf)
         kdetps(jnstru+knstru(jnconf)) = 0
c...  initialize group references
         do 345 iequal=1,nequal
            keqdet((1+maxeq+1)*(iequal-1)+1+maxeq+1) = 0
345      continue
         do 350 icfff=jjcfff,jjcfff-1+kdetps(istruc)
c...  leave out pauli-forbidden determinant
            if (ifindo(krejd,nprejd,idetps(icfff)).gt.0) go to 350
            do 360 iequal=1,nequal
               jieq = (iequal-1)*(1+maxeq+1)
               nieq = keqdet(jieq+1)
               do 370 iieq=1,nieq
                  if (idetps(icfff).eq.keqdet(jieq+1+iieq)) then
c...  current determinant belongs to a group of equal determinants
                     if (keqdet(jieq+1+maxeq+1).eq.0) then
c...  in the current structure the group has not yet been refered to
                        idetps(jncfff+nncfff) = kscr(idetps(icfff))
                        coefff(jncfff+nncfff) = 
     &                  parity(idetps(icfff))*coefff(icfff)
c...  store adress of group on idetps/coefff for current structure
                        keqdet(jieq+1+maxeq+1) = jncfff + nncfff
                        nncfff = nncfff + 1
                        kdetps(jnstru+knstru(jnconf)) = 
     &                  kdetps(jnstru+knstru(jnconf)) + 1
                     else
c...  in the current structure the group has already been refered to
                        coefff(keqdet(jieq+1+maxeq+1)) = 
     &                  coefff(keqdet(jieq+1+maxeq+1)) + 
     &                  parity(idetps(icfff))*coefff(icfff)
                     end if
                     go to 350
                  end if
370            continue
360         continue
c...  current determinant is legal and not in an equals group
            idetps(jncfff+nncfff) = kscr(idetps(icfff))
            coefff(jncfff+nncfff) = parity(idetps(icfff))*coefff(icfff)
            nncfff = nncfff + 1
            kdetps(jnstru+knstru(jnconf)) = 
     &      kdetps(jnstru+knstru(jnconf)) + 1
350      continue
         jjcfff = jjcfff + kdetps(istruc)
         if (kdetps(jnstru+knstru(jnconf)).gt.0) knstru(jnconf) =
     &                                           knstru(jnconf) + 1
340   continue
      call normr(kdetps(jnstru),idetps(jncfff),coefff(jncfff),
     &           knstru(jnconf),nwstru,nwcfff)
c...  number of structures might have been reduced
      if (nwstru.eq.0) then
         null = .true.
         return
      end if
      knstru(jnconf) = nwstru
      nncfff = nwcfff
c...  the excitation process might have generated linear dependent
c...  structures. trace linear dependencies and eliminate redundant
c...  structures. 
      if (knstru(jnconf).gt.1) then
         if (nbeta-ndocc.eq.0) then
            knstru(jnconf) = 1
            nncfff         = 1
         else
c...  calculate mstruc(nsocc) and mdet(nsocc) (conform initsy)
            spin = (mult-1.0d0)/2.0d0
            n2s  = (nsocc/2.0d0-spin+.0001d0)
            if (n2s.eq.0) then
               mstruc = 1
            else
               new1   = newtont(nsocc,n2s)
               n2s1   = n2s-1
               new2   = newtont(nsocc,n2s1)
               mstruc = new1 - new2
            end if
            if(mstruc.gt.maxspvb)call vberr('maxspvb exceeded (excits)')
            nalph = (nsocc + mult - 1)/2
            mdet  = newtont(nsocc,nalph)
c...  dynamical core allocation for lindep 
            ksocc  = 1         
            kdetp  = ksocc  +   nsocc
            idetp  = kdetp  +   mstruc
            kkkdet = idetp  +   mstruc*mdet
            kkdete = kkkdet +   mdet*nelec
            idetnr = kkdete +   mdet*nelec
            ksbeta = idetnr +   mdet
            kendin = ksbeta +  nsbeta          
c      
            kcoeff = 1
            kcoec  = kcoeff + mstruc*mdet 
            kendre = kcoec  + mstruc*mdet
c
            if (kendin.gt.mscr) call vberr 
     &         ('core fault in doubst before call to lindep')
c...  collect single occupieds and determinants for lindep
            call icopy(nalpha-ndocc,kdeter(1,1),1,kscr(1),1)
            call icopy(nbeta -ndocc,kdeter(nalpha+1,1),1,
     &                 kscr(nalpha-ndocc+1),1)
            call ordernumbs(kscr(ksocc),nsocc)
            call icopy(kntdet(jnconf)*nelec,kdeter,1,kscr(kkdete),1)
            call lindep(jnconf,1,kscr(ksocc),nsocc,mdet,mstruc,
     &                  jnstru,jncfff,dummy,knstru,kntdet,kdetps,
     &                  idetps,coefff,dummy,newstr,newdet,newcff,
     &                  kscr(kdetp),kscr(idetp),scr(kcoeff),
     &                  kscr(kkkdet),melec1,.false.,kscr(kkdete),
     &                  kscr(idetnr),scr(kcoec),kscr(ksbeta))
            knstru(jnconf) = newstr
            kntdet(jnconf) = newdet
            nncfff = newcff
            call icopy(newstr,kscr(kdetp) ,1,kdetps(jnstru),1)
            call icopy(newcff,kscr(idetp) ,1,idetps(jncfff),1)
            call fmove(scr(kcoeff),coefff(jncfff),newcff)
            call icopy(newdet*nelec,kscr(kkkdet),1,kdeter        ,1)
         end if
      end if
c
      return
c
4     continue
c...
c...  jcase 4: virtual is singly occupied, active is doubly occupied
c...
c
c...
c...  copy and modify determinants
c...
c
      jposa = nalpha - ndocc + (jposa-nsocc+1)/2
c...  use educated guess for jposv at first determinant, based on
c...  ordering determinants after spin projection
      if (jposv.le.nbeta-ndocc) then
         jposv = nalpha + jposv
      else
         jposv = jposv - (nbeta-ndocc)
      end if
      do 400 idet=1,kntdet(jrconf)
c... search for position (jposv) of jvirt on column idet of kdeter
c... guess for jposv if idet exceeds 1: jposv(idet) = jposv(idet-1)
         call ifindd(kdeter(1,idet),jvirt,jposv,nalpha,ndocc)
c...  excite determinant
         nposv = nalpha - ndocc + irank(kdeter(nalpha-ndocc+1,
     &                                  idet),ndocc,jvirt,
     &                                  jposa-(nalpha-ndocc),0)
c...  mind interfering permutations
         if (jposv.le.nalpha) then
            nposa =          irank(kdeter(       1,idet),
     &                             nalpha-ndocc,jactiv,jposv,0)
            nperm = iabs((nposv+nbeta)-(jposa+nbeta)) +
     &              (iabs(nposa-jposa) + 
     &               iabs(nposv-jposv) - 1) 
         else
            nposa = nalpha + irank(kdeter(nalpha+1,idet),
     &                             nbeta -ndocc,jactiv,jposv,0)
            nperm = iabs(nposv-jposa) +
     &              (iabs(nposa-(jposa+nbeta)) + 
     &               iabs((nposv+nbeta)-jposv) - 1)
         end if
         call open1(kdeter(1,idet),nposa,jposv)
         kdeter(nposa,idet) = jactiv
         call open1(kdeter(1,idet),nposv      ,jposa)
         kdeter(nposv      ,idet) = jvirt
         call open1(kdeter(1,idet),nposv+nbeta,jposa+nbeta)
         kdeter(nposv+nbeta,idet) = jvirt
         parity(idet) = (-1.0d0)**nperm
400   continue
      kntdet(jnconf) = kntdet(jrconf)
c...
c...  copy structure-definitions
c...
      call icopy(nrcfff,idetps(jrcfff),1,idetps(jncfff),1)
      do 410 icfff=1,nrcfff
         coefff(jncfff-1+icfff) = 
     &   parity(idetps(jncfff-1+icfff))*coefff(jrcfff-1+icfff)
410   continue
      nncfff  = nrcfff
      call icopy(knstru(jrconf),kdetps(jrstru),1,kdetps(jnstru),1)
      knstru(jnconf) = knstru(jrconf)
c
      return
c
c...
c...        ------------------------------------------------------------
c...
c
      end
      integer function ialign(iref,ich,itake,nd)
      implicit REAL  (a-h,o-z) ,  integer   (i-n)
c
c...   align ich so that it is identical to iref and take itake along
c...   ialign keeps the sign of the total permutation
c
      dimension iref(nd),ich(nd),itake(nd)
c
      iperm = 1
      do 100 ir=1,nd-1
      if (iref(ir).ne.ich(ir)) then
c...      look for iref
         do 10 ic=ir+1,nd
10       if (ich(ic).eq.iref(ir)) go to 20
         call vberr(' impossible task for ialign')
c....    interchange ic and ir
20       iperm = -iperm
         isave = ich(ir)
         ich(ir) = ich(ic)
         ich(ic) = isave
         isave = itake(ir)
         itake(ir) = itake(ic)
         itake(ic) = isave
      end if
100   continue
c
      ialign = iperm
c
      return
      end
      integer function ifind(knumbs,nnumbs,jtest)
c...
c... return position of jtest on knumbs
c...        0 if jtest is not present
c... no assumptions on the structure of knumbs are made
c... 
      implicit REAL (a-h,o-z), integer (i-n)
c
      dimension knumbs(*)
c
      do 10 i=1,nnumbs
         if (jtest.eq.knumbs(i)) then
            ifind = i 
            return
         end if
10    continue
      ifind = 0
c
      return
      end
      subroutine ifindd(kdeter,jorb,jpos,nalpha,ndocc)
c...  completely devoted to excits
c...
c...  find position of jorb on kdeter (jpos). jorb must be in 
c...  in the singles-part of eighter the alpha or the beta-block.   
c...  the value of jpos at input is used as guess, so it must be part
c...  of eighter one of the blocks.
c...
c...  assumed construction of determinants:
c...       /        alpha block        /         beta block        /
c...       /   singles   /   doubles   /   singles   /   doubles   /
c...  size: nalpha-ndocc     ndocc      nbeta-ndocc      ndocc
c...
c...  within each block, orb(i+1) >= orb(i) holds. 
c...
c
      implicit REAL (a-h,o-z), integer (i-n)
c
      common /indep/ time0,nelec,mult
c
      dimension kdeter(*)
c
      if (jpos-(nalpha-ndocc)) 10,10,20
c...  alpha-block 
10    if (kdeter(jpos).eq.jorb) return
      if (kdeter(jpos).gt.jorb) then
         if (jpos-1.eq.0             .or.kdeter(jpos-1).lt.jorb) then
            jpos = nalpha + 1
            go to 20
         else
            jpos = jpos - 1
         end if
      else
         if (jpos+1.eq.nalpha-ndocc+1.or.kdeter(jpos+1).gt.jorb) then
            jpos = nalpha + 1
            go to 20
         else
            jpos = jpos + 1
         end if
      end if
      go to 10     
c... beta-block 
20    if (kdeter(jpos).eq.jorb) return
      if (kdeter(jpos).gt.jorb) then
         if (jpos-1.eq.nalpha        .or.kdeter(jpos-1).lt.jorb) then
            jpos = 1
            go to 10
         else
            jpos = jpos - 1
         end if
      else
         if (jpos+1.eq.nelec-ndocc+1 .or.kdeter(jpos+1).gt.jorb) then
            jpos = 1
            go to 10
         else                         
            jpos = jpos + 1
         end if
      end if
      go to 20
c
      end
      integer function ifindo(knumbs,nnumbs,jtest)
c...
c...  return position of jtest on ordered array knumbs
c...         0 if jtest is not present
c...  it is asumed that knumbs(i+1) > knumbs(i) for all i
c...
      implicit REAL (a-h,o-z), integer (i-n)
c
      dimension knumbs(*)
c
      if (nnumbs.le.0) then
         ifindo = 0
         return
      end if
c
      ifindo = 1 
10    if (knumbs(ifindo).eq.jtest) return
      if (knumbs(ifindo).lt.jtest.and.ifindo+1.le.nnumbs) then
         ifindo = ifindo + 1
         go to 10
      else
         ifindo = 0
      end if 
c 
      return
      end
      subroutine initsy(kconf,melec1,nconf,nspath,msp)
c...
c...  'initsy' puts the factorials 0 to mmsocc at common block /fac/,
c...  and determines vital maxima for partitioning symbolic
c...  input : msocc (/hans/) = max. singly occupieds ; mult = 
c...                                              multiplicity (/indep/)
c...  ** called by cresti **
c...
      implicit REAL  (a-h,o-z) , integer (i-n)
c
INCLUDE(common/turtleparam)
INCLUDE(../m4/common/iofile)
c
      common /indep/  time0,nelec,mult
      common /fac/    kfacul(0:mmsocc)
      common /hans/   msocc,mstrco,mdetco,mdetst,mperm,mstruc
INCLUDE(common/logics)
c
      dimension kconf (melec1,*),
     &          nspath(*)
c
c...  if iijj... --> ikjl... nsocc = nsocc + 4
      if (mrsd) msocc = min(msocc+4,nelec)    
c
      if (msocc.gt.mmsocc) call dimerr('single occupied orbitals',
     &                                 mmsocc)
c...
c...                       calculate factorials
c...
      kfacul(0) = 1
      do 10 ifac=1,mmsocc
         kfacul(ifac) = kfacul(ifac-1) * ifac
10    continue
c
c...  calculate  spin-possibilities and  determinants
c...  for multiplicity mult and msocc singly occupieds
c
     
      spin = (mult-1.0d0)/2.0d0
c...   n2s is  linked pairs
      n2s = (msocc/2.0d0-spin+.0001d0)
c...   mstrco =  spin-possibilities =  structures/config
      if (n2s.eq.0) then
         mstrco = 1
      else
c
         if (msp.eq.0) then
            new1 = newtont(msocc,n2s)
            n2s1 = n2s-1
            new2 = newtont(msocc,n2s1)
            mstrco = new1 - new2
            msp = mstrco
         else
            mstrco = msp
         endif
         if (msp.gt.maxspvb) then
            write(iwr,600) msp,maxspvb
600         format(
     o         /'   **************************************************',
     1         /'   ** max. # spin functions possible ',i10,   '    **',
     2         /'   ** max. # spin functions allowed  ',i10,   '    **',
     3         /'   ** possibly adjust your problem or turtleparam  **',
     4         /'   **************************************************',
     5         /' ')
            msp = maxspvb
            mstrco = msp
         end if
c
      end if
c
      nalpha = (msocc + mult - 1)/2
c...   mdetst =  determinants / rumer structure
      mdetst = 2**((msocc-mult+1)/2)
c
c...   mdetco =  determinants / configuration
c
      if (msp.eq.0) then
         mdetco = newtont(msocc,nalpha)
      else
         mdetco = msp*mdetst
      endif
c
c...    mperm = n2s over (n2s/2)
c
      n2s2 = n2s / 2
      mperm = newtont(n2s,n2s2)
c
c...  total number of structures for known configurations
      mstruc = 0
      do 20 iconf=1,nconf
         if (nspath(iconf).eq.0) then
c...  calculate number of structures (all structures are used)
            ndocc = kconf(nelec+1,iconf)
            nsocc = nelec - 2*ndocc
c...  npair is  linked pairs
            npair = (nsocc/2.0d0-spin+.0001d0)
            if (npair.eq.0) then
               nstrco = 1
            else
               new1 = newtont(nsocc,npair)
               n2s1 = npair-1
               new2 = newtont(nsocc,n2s1)
               nstrco = new1 - new2
            end if
         else
            nstrco = nspath(iconf)
         end if
         nstrco = min(nstrco,mstrco)
         mstruc = mstruc + nstrco
20    continue
c...  check if only rumers are wanted
      allrum = .true.
      do 30 iconf=1,nconf
         if (.not.rumer(iconf)) then
            allrum = .false.
            go to 40
         end if
30    continue
40    continue
c
      return
      end
      integer function irank(knumbs,nnumbs,jtest,neglc1,neglc2)
c...
c...  return the position where jtest can be put so that it's 
c...  presence will not disturb the ordering of knumbs. (assume a right
c...  shift is performed before excitation)
c...  the ordering is characterized by: knumbs(i+1) > knumbs(i)
c...  neglect the contents of elements knumbs(neglc1) and knumbs(neglc2)
c...
      implicit REAL (a-h,o-z), integer (i-n)
c
      dimension knumbs (*)
c
      irank = 1
      do 10 itest=1,nnumbs
         if (itest.eq.neglc1.or.itest.eq.neglc2) then 
            continue 
         else if (knumbs(itest).le.jtest) then
            irank = irank + 1
         else
            return
         end if 
10    continue
c 
      return
      end
      function isconf (karray,nelem)
c***********************************************************************
c
      implicit REAL  (a-h,o-z) , integer (i-n)
c
c***********************************************************************
c
c    'isconf' determines the symmetry-number of the configuration that
c     is defined by     'karray'.
c
      dimension karray ( * )
c
INCLUDE(common/turtleparam)
      common/vbtable/nirr,mult(8,8),isymr(2,8),jgroup,ifrep,iro(mxorbvb)
c
      ii = iro(karray(1))
      do 10 i=2,nelem
         ii = mult(ii,iro(karray(i)))
10    continue
      isconf = ii
c
      return
      end
      subroutine lindep(keqcf ,neqcf ,ksocc ,nsocc ,mdet,mstruc,
     &                  jstruc,jcfff ,jpacd ,
     &                  knstru,kntdet,kdetps,idetps,coefff,pacdet,
     &                  nnstru,nndet ,nncfff,kdetp ,idetp ,coeff ,
     &                  kdeter,melec1,pac   ,
     &                  kkdete,idetnr,coec,ksbeta)
c...
c...  'lindep' selects a set of linear independent structures from a 
c...  group of equal configurations.
c...  linear dependencies are detected by means of schmidt-
c...  orthogonalization.
c...
c...  the number of single occupied beta orbitals must exceed 0!
c...
c...  definitions of selected set: nnstru, nndet, nncfff
c...                               kdetp, idetp, coeff, kdeter
c...
      implicit REAL (a-h,o-z), integer (i-n)
c
INCLUDE(common/c8_16vb)
c
      logical pac
c
      common /indep/ time0,nelec,mult
c
      dimension keqcf (*),
c
     &          knstru(*),
     &          kntdet(*),
     &          kdetps(*),
     &          idetps(*),
     &          coefff(*),
     &          pacdet(*),
c
     &          kdetp (*),
     &          idetp (*),
     &          coeff (*),
     &          kdeter(melec1-1,*),

     &          ksocc (*),
     &          kkdete(melec1-1,*),
     &          idetnr(*),
     &          coec  (mdet,*),
     &          ksbeta(*)
INCLUDE(common/vbcri)
c      
      nalpha = (nelec + mult - 1)/2
      nbeta  = nelec - nalpha
      ndocc  = (nelec-nsocc)/2
      nsbeta = nbeta - ndocc    
c...  set null-label for resultant determinants
      do 10 idet=1,mdet
         kdeter(1,idet) = 0
10    continue
c...  initialize size parameters
      nnstru = 0
      nncfff = 0
c...  loop over configuration group
      do 100 ieqcf=1,neqcf     
c...
c...  set definitions for current configuration
c...
         if (ieqcf.gt.1) then
c...  determine start adresses and unpack determinants
            do 110 iconf=keqcf(ieqcf-1),keqcf(ieqcf)-1
               do 120 istruc=jstruc,jstruc-1+knstru(iconf)
                  jcfff = jcfff + kdetps(istruc)
120            continue 
               jstruc = jstruc + knstru(iconf)
               if (pac) jpacd = jpacd + 
     &                  (nelec*kntdet(iconf)-1)/(64/n8_16) + 1
110         continue
         end if
         if (pac) 
     &   call unpack(pacdet(jpacd),n8_16,kkdete,
     &               nelec*kntdet(keqcf(ieqcf)))
         do 130 idet=1,kntdet(keqcf(ieqcf))
c...  gather singly occupied beta's
            jsbeta = 1
            do 140 ibeta=nalpha+1,nelec
               if (ifind(ksocc,nsocc,kkdete(ibeta,idet)).gt.0) then
                  ksbeta(jsbeta) = kkdete(ibeta,idet)
                  jsbeta = jsbeta + 1
               end if
140         continue
            call ordernumbs(ksbeta,nsbeta) 
c...  lay connection between current ordering determinants and their 
c...  logical number, based upon single occupied beta orbitals
            idetnr(idet) = lognr(ksbeta,nsbeta,1,ksocc,nsocc)
130      continue      
c...  initialize coefficients       
         do 170 istruc=nnstru+1,nnstru+knstru(keqcf(ieqcf))
            do 180 idet=1,mdet
               coec(idet,istruc) = 0.0d0
180         continue
170      continue      
c...
c...  schmidt-orthogonalize structures one by one to (survived) 
c...  previously schmidt-orthogonalized structures. reject structure
c...  if orthogonalization turns it into a null-vector. store original 
c...  definitions of structures that have proven to be relevant.
c...
         jjcfff = 0
         do 200 istruc=jstruc,jstruc-1+knstru(keqcf(ieqcf))
c...  redefine structure in terms of logically ordered determinants
            do 210 icfff=jcfff+jjcfff,jcfff-1+jjcfff+kdetps(istruc)
               coec(idetnr(idetps(icfff)),nnstru+1) = coefff(icfff)
210         continue
c...  orthogonalize
            do 220 iprevs=1,nnstru
               prod = 0.0d0
               do 230 idet=1,mdet
                 prod = prod + coec(idet,nnstru+1)*coec(idet,iprevs)
230            continue
               do 240 idet=1,mdet
                  coec(idet,nnstru+1) = coec(idet,nnstru+1) - 
     &                                  prod*coec(idet,iprevs)
240            continue
220         continue
c...  normalize
            cnorm = 0.0d0
            do 250 idet=1,mdet
               cnorm = cnorm + coec(idet,nnstru+1)**2
250         continue
            if (cnorm.gt.cridep) then
c...  current structure is a relevant extension to the previously  
c...  defined set of structures
               cnorm = 1.0d0/dsqrt(cnorm)
               do 260 idet=1,mdet
                  coec(idet,nnstru+1) = coec(idet,nnstru+1)*cnorm
260            continue
c...  store non-orthogonalized version of current structure 
               do 270 icfff=1,kdetps(istruc)
                  idetp(nncfff+icfff) = 
     &                              idetnr(idetps(jcfff-1+jjcfff+icfff))
                  coeff(nncfff+icfff) = coefff(jcfff-1+jjcfff+icfff)
270            continue
               kdetp(nnstru+1) = kdetps(istruc)
c...  check if any new determinants must be defined
               do 280 icfff=1,kdetps(istruc)
                  if (kdeter(1,idetp(nncfff+icfff)).eq.0) then
c...  determinant has not yet been defined. do it now
                     call icopy(nelec,
     &                          kkdete(1,idetps(jcfff+jjcfff-1+icfff)),
     &                          1,kdeter(1,idetp(nncfff+icfff)),1)
                  end if
280            continue
               nncfff = nncfff + kdetp(nnstru+1)
               nnstru = nnstru + 1
               if (nnstru.eq.mstruc) then
c...  a complete set of spin functions has been defined
                  nndet = mdet
                  return
               end if
            end if
            jjcfff = jjcfff + kdetps(istruc)
200      continue
100   continue 
c...
c...  remove redundant determinants
c...
      nrejd = 0
      do 300 idet=1,mdet
         if (kdeter(1,idet).eq.0) then
            do 310 icfff=1,nncfff
               if (idetp(icfff).gt.idet-nrejd) idetp(icfff) = 
     &                                         idetp(icfff) - 1
310         continue
            nrejd = nrejd + 1
         else
            call copi(kdeter(1,idet),kdeter(1,idet-nrejd),nelec,1)
         end if
300   continue
      nndet = mdet - nrejd
c
      return
      end
      integer function lognr(karray,npos,incr,knumbs,nnumbs)
c...
c...  return logical number to a row of numbers
c...
c...  the row is defined by karray(1+i*incr), i=1,npos.
c...  karray(1+(i+1)*incr) > karray(1+i*incr) for every i.
c...  each element of karray is member of knumbs. for every i,
c...  knumbs(i+1) > knumbs(i) holds . consequently, the position of
c...  an element of karray on knumbs gives the rank of this number among
c...  the set of numbers held by knumbs. it is the ensemble of ranks of
c...  elements of karray that the logical number is based on.
c...  be converted that the logical number is based on.
c...
c...  suppose there are 5 positions and 7 numbers. then the order in 
c...  which all possible combinations can be realized systematically is
c...  as follows: (numbers refering to ranks of elements)
c...                         1 2 3 4 5
c...                         1 2 3 4 6 
c...                         1 2 3 4 7
c...                         1 2 3 5 6
c...                         1 2 3 5 7
c...                         1 2 3 6 7
c...                         . . . . .
c...                         3 4 5 6 7
c...  the first row would be given logical number 1, etc..
c...
c...  now consider how the logical number of a general row, consisting
c...  of npos numbers who are member of a set of nnumbs numbers, is
c...  build up. the ranks of the elements of the row are considered 
c...  piece by piece. the general filosophy is that on going from one 
c...  element to another some possibilities may be skipped and we can
c...  raise the minimal logical number. this goes on and on till we have
c...  reached the last element. 
c...
c...  suppose the rank of the first element exceeds 1. then the
c...  possibilities that are skipped are:
c...
c...  - rank karray(1) : 1
c...    ranks of karray(2) to karray(npos) : all combinations of numbers
c...                                         2 to nnumbs 
c...    = nnumbs-1 over npos-1 possibilities
c...  - rank karray(1) : 2
c...    ranks of karray(2) to karray(npos) : all combinations of numbers
c...                                         3 to nnumbs 
c...    = nnumbs-2 over npos-1 possibilities
c...  - ...
c...  - rank karray(1) : karray(1)-1
c...    ranks of karray(2) to karray(npos) : all combinations of numbers
c...                                         karray(1) to nnumbs
c...    = nnumbs-karray(1) over npos-1 possibilities
c...
c...  the number of apperently skipped combinations after revealing the
c...  rank of the first element is thus calculated. if we go on to the 
c...  second element than we are dealing with skipped combinations if
c...  irank karray(2) > irank karray(1)+1. the number of skipped
c...  combinations is calculated in simular way to what has been shown 
c...  above. this goes on and on till the last element has been reached
c...
c...  if we define rank(karray(0)) = 0, then we can develop the
c...  following expression:
c...  lognr = 1 + sigma(ipos=1,npos) 
c...              sigma(j=rank(karray(ipos-1))+1,rank(karray(ipos))-1
c...              nnumbs-j over npos-ipos
c...
      implicit REAL (a-h,o-z), integer (i-n)
c
      dimension karray(*),
     &          knumbs(*)
c
      if (npos.eq.0) then
         lognr = 0
         return
      end if
c
      lognr  = 1
      ipos   = 1
      joccur = ifindo(knumbs,nnumbs,karray(ipos*incr))
      do 10 j=0+1,joccur-1
         lognr = lognr + newtont(nnumbs-j,npos-ipos) 
10       continue
         do 20 ipos=2,npos
            joccm1 = joccur
            joccur = ifindo(knumbs,nnumbs,karray(ipos*incr))
            do 30 j=joccm1+1,joccur-1
               lognr = lognr + newtont(nnumbs-j,npos-ipos) 
30          continue
20       continue
c
      return
      end
      subroutine maxi_crestr(korb,ninner,kvirt,nvirt,kntdet,knstru,
     &                       nconf,imax,maxdet,maxstr)
c...
c...  calculate some maxima
c...
c...  originally subr max, max is intrinsic function on SiliconGraphics
c...  maxi gave clash with peigss
c...
      implicit REAL (a-h,o-z), integer (i-n)
c
      common /indep/ time0,nelec,mult
c
      dimension korb  (*),
     &          kvirt (*),
     &          kntdet(*),
     &          knstru(*) 
c
c...  highest orbital number
      imax = korb(ninner)
      do 10 ivirt=1,nvirt
         imax = max(imax,kvirt(ivirt))
10    continue
c...  maximum number of determinants/structures per configuration
      maxdet = kntdet(1)
      maxstr = knstru(1)
      do 20 iconf=2,nconf
         maxdet = max(maxdet,kntdet(iconf))
         maxstr = max(maxstr,knstru(iconf))
20    continue
c
      return
      end
      subroutine mrsdd(kconf ,melec1,nsp0  ,nspath,ispath,maxsp , 
     &                 korb  ,ninner,kvirt ,nvirt ,idconf,
     &                 knstru,kntdet,kdetps,idetps,coefff,pacdet,
     &                 kstruc,coeff ,kdeter,klogic,parity,kscr  ,
     &                 mscr  ,scr   ,mconf ,irumer)
c...
c...  driver for multi-reference singles and doubles    
c...
c...  eighter one of two types of excitations is performed (depending on
c...  logical variable excstr) :
c...  - on configuration basis: a spin projection is performed after the
c...    excitation has taken place. selection of spin functions is not
c...    allowed. 
c...  - on structure basis: the excited structures are builded directly 
c...    by means of a spin-adapted excitation operator. selection of 
c...    spin functions may take place at defining the reference space.
c...
c...  external virtuals are defined symbolically. orbitals 254 and 255
c...  represent them. now really 2**n8_16-2 and 2**n8_16-1
c...
c...  the performance of the excitation processes is increased by 
c...  thoroughly presuming the ordering of spin blocks.
c...
c...  sequencially performed tasks:
c...  - spin projection of reference configurations 
c...  - internal ordering of reference configurations/determinants
c...  - definition of single excited structures
c...  - definition of double excited structures
c...  - removal of redundants
c...  - construction of data tape
c...
      implicit REAL (a-h,o-z), integer (i-n)
c
INCLUDE(common/c8_16vb)
INCLUDE(common/turtleparam)
c
      dimension kconf  (melec1,*)
      dimension nspath (*),
     &          ispath (maxsp,*),    
c
     &          korb   (*),
     &          kvirt  (*),
     &          idconf (4,*),
c
     &          knstru (*),
     &          kntdet (*),
     &          kdetps (*),
     &          idetps (*),
     &          coefff (*),
     &          pacdet (*),
c
     &          kstruc (*),
     &          coeff  (*),
     &          kdeter (melec1-1,*),
     &          klogic (*),
     &          parity (*),
     &          kscr   (*),
     &          irumer (*),
     &          scr    (*),
     &          idoubz(1)
c
INCLUDE(common/logics)
      common /indep/  time0,nelec,mult
      common /hans/   msocc,mstrco,mdetco,mdetst,mperm,mstruc
      common /size/   nrefco,nsinin,nsinex,ndouin,ndoumi,ndoue1,ndoue2
INCLUDE(../m4/common/iofile)
INCLUDE(../m4/common/timeperiods)
c
      logical nofail
c                            
      write(iwr,*) 'input-arrays in mrsdd'
      write(iwr,*) 'kconf  ', ((kconf(j,i) , j=1,nelec), i=1,nrefco)   
      write(iwr,*) 'nspath ',  (nspath(i),   i=1,nrefco)
      write(iwr,*) 'ispath ',  (ispath(1,i), i=1,nrefco)
      write(iwr,*) 'kvirt  ',  (kvirt(i),    i=1,nvirt)
      write(iwr,*) 'korb  ',   (korb(i),     i=1,ninner+2) 
      write(iwr,*) 'idconf ', ((idconf(j,i), j=1,4), i=1,nrefco)    
      write(iwr,*) 'rumer  ',  (rumer(i),    i=1,nrefco)
      write(iwr,5)
5     format(///1x,70(1h*),/,
     &         13x,'------ generation of singles and doubles ------')
      write(iwr,15)
15    format(/ 13x,',       (symbolic definition of externals)')
      if (excstr) then
         write(iwr,25) 'structure'
      else
         write(iwr,25) 'configuration'
      end if   
25    format(/ 13x,'       - excitation on ',a,' basis')
      write(iwr,35) nrefco
35    format(//13x,'               # reference configurations: ',i4,/,
     &       1x,70(1h*))
c
      nalpha = (nelec + mult - 1)/2
      nbeta  = nelec - nalpha
      norb   = ninner + 2
      if (excstr) then
         nofail = .true.   
      else
         nofail = .false.
      end if
c...
c...  perform spinprojection of reference configurations           
c...
      krejcf = 1
      kkscr  = krejcf + nrefco 
c...  kendin = kkscr  + mspine
c...  if (kendin.gt.mscr) call vberr(' insufficient core memory at call t
c... &o projc')
      call projc(kconf,melec1,nrefco,1,nspath,ispath,maxsp,nofail,
     &           knstru,kntdet,kdetps,idetps,coefff,pacdet,nstruc,ncfff,
     &           npacd,ndet,kscr(krejcf),nrejcf,kstruc,coeff,kdeter,
     &           klogic,kscr(kkscr),mscr,scr,irumer)
c...  remove configurations leading to nul-functions from idconf 
      if (nrejcf.gt.0) call skipi(idconf,4,nrefco,kscr(krejcf),nrejcf,
     &                            idconf)
      nconf = nrefco - nrejcf
      time1 = cpulft(1)
c...
c...  complete ordering of spin blocks (they are ordered already for
c...  configurations without structure selection) 
c...
      if (excstr) then
c...  order reference configurations for excitation on structure basis
         do 10 iconf=nsp0+1,nrefco
            nsocc = nelec - 2*kconf(nelec+1,iconf)
            call orders(kconf(1,iconf),melec1,1,1,nsocc)
10       continue
c...  order spin blocks of reference determinants for excitation on
c...  structure basis   
         jpacd  = 1
         jstruc = 1
         jcfff  = 1
         do 20 iconf=1,nsp0
            jpacd = jpacd + (nelec*kntdet(iconf)-1)/(64/n8_16) + 1
            do 30 istruc=jstruc,jstruc-1+knstru(iconf)
               jcfff = jcfff + kdetps(istruc)
30          continue
            jstruc = jstruc + knstru(iconf)
20       continue
         do 40 iconf=nsp0+1,nrefco
            ndocc = kconf(nelec+1,iconf)
            call unpack(pacdet(jpacd),n8_16,kdeter,nelec*kntdet(iconf))
            do 50 idet=1,kntdet(iconf)
               call orderp(kdeter(1       ,idet),nalpha-ndocc,nperm1)
               call orderp(kdeter(nalpha+1,idet),nbeta -ndocc,nperm2)
               parity(idet) = (-1.0d0)**(nperm1+nperm2)
50          continue
            call pack(pacdet(jpacd),n8_16,kdeter,nelec*kntdet(iconf))
            jpacd = jpacd + (nelec*kntdet(iconf)-1)/(64/n8_16) + 1
c...  correct coefficients for parity changes due to permutations
            do 60 istruc=jstruc,jstruc-1+knstru(iconf)
               do 70 icfff=jcfff,jcfff-1+kdetps(istruc)
                  coefff(icfff) = parity(idetps(icfff))*coefff(icfff)
70             continue
               jcfff = jcfff + kdetps(istruc)
60          continue
            jstruc = jstruc + knstru(iconf)
40       continue
      end if
      call pri(kconf,melec1,idconf,knstru,kntdet,kdetps,idetps,
     &         coefff,pacdet,kdeter,.false.,nconf,nconf)
      if (excstr) then
c...  partition kscr for sings/doubs
         kcurcf = 1
         krejd  = kcurcf + nelec + 1
         keqdet = krejd  + mdetco
         kkkscr = keqdet + mdetco*(mdetco+2)
      else
c...  partition kscr for singc/doubc
         kcurcf = 1
         kkscr  = kcurcf + nelec + 1
      end if
c...
c...  ++++++++++++++++++++++++++ s i n g l e s +++++++++++++++++++++++++
c...
      do 100 ivirt=1,ninner+1
         if (ivirt.eq.ninner+1) nsinin = nconf - nrefco
         do 110 iact=1,ninner
            if (iact.eq.ivirt) go to 110
            jvirt  = korb(ivirt)
            jactiv = korb(iact)     
cccc            write(iwr,*) jactiv, ' to ', jvirt
            if (excstr) then
c...  excite on structure level
               call sings(kconf,melec1,nrefco,korb,norb,idconf,
     &                    jactiv,jvirt,knstru,kntdet,kdetps,idetps,
     &                    coefff,pacdet,nconf,nstruc,ncfff,npacd,ndet,
     &                    kdeter,parity,kscr(krejd),kscr(keqdet),
     &                    kscr(kkkscr),scr,kscr(kcurcf),mscr-kkkscr+1,
     &                    mconf)
            else
c...  excite on configuration level  
               call singc(kconf,melec1,nrefco,nspath,ispath,
     &                    maxsp,korb,norb,idconf,jactiv,jvirt,knstru,
     &                    kntdet,kdetps,idetps,coefff,pacdet,nconf,
     &                    nstruc,ncfff,npacd,ndet,kstruc,coeff,kdeter,
     &                    klogic,parity,kscr(kkscr),kscr(kcurcf),
     &                    mscr-kkscr+1,scr,mconf,irumer)
            end if
110      continue
100   continue
      nsinex = nconf - nsinin - nrefco
c...
c...  ++++++++++++++++++++++++++ d o u b l e s +++++++++++++++++++++++++
c...
      nprevc = nconf
      do 200 itype=1,6
         if (itype.eq.1) then
c...  internals/internals
            jf1 = 1
            jl1 = ninner
            jf2 = 1
            jl2 = ninner
         else if (itype.eq.2) then
            ndouin = nconf - nprevc
            nprevc = nconf
c...  internals/externals
            jf1 = 1
            jl1 = ninner
            jf2 = ninner + 1
            jl2 = ninner + 1
         else if (itype.eq.3) then
c...  externals/internals
            jf1 = ninner + 1
            jl1 = ninner + 1
            jf2 = 1
            jl2 = ninner
         else if (itype.eq.4) then
            ndoumi = nconf - nprevc
            nprevc = nconf
c...  externals/externals (not equal)
            jf1 = ninner + 1
            jl1 = ninner + 1
            jf2 = ninner + 2
            jl2 = ninner + 2
         else if (itype.eq.5) then
c...  externals/externals (not equal) (=type4 reversed)
            jf1 = ninner + 2
            jl1 = ninner + 2
            jf2 = ninner + 1
            jl2 = ninner + 1
         else
            ndoue1 = nconf - nprevc
            nprevc = nconf
c...  externals/externals (equal)
            jf1 = ninner + 1
            jl1 = ninner + 1
            jf2 = ninner + 1
            jl2 = ninner + 1
         end if
         do 210 ivirt1=jf1,jl1
            do 220 ivirt2=jf2,jl2
               do 230 iact1=1,ninner
                  if (iact1.eq.ivirt1) go to 230
                  do 240 iact2=iact1,ninner
                     if (iact2.eq.ivirt2) go to 240
                     jact1  = korb(iact1)
                     jvirt1 = korb(ivirt1)
                     jact2  = korb(iact2)
                     jvirt2 = korb(ivirt2)      
ccccc                     write(iwr,*) jact1, ' to ', jvirt1, ' and '
ccccc                     write(iwr,*) jact2, ' to ', jvirt2   
                     if (excstr) then
c...  excite on structure level
                        call doubs(kconf,melec1,nrefco,korb,norb,idconf,
     &                             jact1,jvirt1,jact2,jvirt2,knstru,
     &                             kntdet,kdetps,idetps,coefff,pacdet,
     &                             nconf,nstruc,ncfff,npacd,ndet,kdeter,
     &                             parity,kscr(krejd),kscr(keqdet),
     &                             kscr(kkkscr),scr,kscr(kcurcf),
     &                             mscr-kkscr+1,mconf)
                     else
c...  excite on configuration level
c...  skip redundant excitations (based upon completeness spin space)
                        if ((ivirt2.ge.ivirt1).and.
     &                      (iact2 .ne.ivirt1).and.
     &                      (ivirt2.ne.iact1 ))     then
                           call doubc(kconf,melec1,nrefco,nspath,
     &                                ispath,maxsp,korb,norb,idconf,
     &                                jact1,jvirt1,jact2,jvirt2,knstru,
     &                                kntdet,kdetps,idetps,coefff,
     &                                pacdet,nconf,nstruc,ncfff,npacd,
     &                                ndet,kstruc,coeff,kdeter,klogic,
     &                                parity,kscr(kkscr),kscr(kcurcf),
     &                                mscr-kkscr+1,scr,mconf,irumer)
                        end if
                     end if
240               continue
230            continue
220         continue
210      continue
200   continue
      ndoue2 = nconf - nprevc
      time2 = cpulft(1)
      write(iwr,45) cpulft(1) - time1
45    format(/,t10,'** configuration generation: ',f6.2,' cpu-seconds', 
     &                                                           ' **')
c...
c...  remove redundant structures/configurations   
c...
      write(iwr,*)'nconf,nrefco,nsinin,nsinex,ndouin,ndoumi,
     &             ndoue1,ndoue2'
      write(iwr,*) nconf,nrefco,nsinin,nsinex,ndouin,ndoumi,
     &             ndoue1,ndoue2
      krejcf = 1
      kkscr  = krejcf + nconf
      if (kkscr.gt.mscr) call vberr(' core fault before call to remred')
      call remred(idconf,knstru,kntdet,kdetps,idetps,coefff,pacdet,
     &            nconf,nstruc,ncfff,npacd,ndet,kscr(krejcf),
     &            kscr(kkscr),scr,melec1,mscr-kkscr+1)        
c...  calculate explicit number of configurations per type-for output 
      n1 = nsinin
      n2 = nsinex*nvirt
      n3 = ndouin
      n4 = ndoumi*nvirt
      n5 = ndoue1*(nvirt)*(nvirt-1)/2 + ndoue2*nvirt    
      write(iwr,55) cpulft(1) - time2
55    format(/,t10,'** removal of redundants:    ',f6.2,' cpu-seconds',
     &                                                           ' **')   
      write(iwr,65) nconf,nstruc,ndet
65    format(//8x,57(1h*),/,t30,'final results',//,
     &          t24,'number of configurations : ',i6,/
     &          t24,'number of structures     : ',i6,/
     &          t24,'number of determinants   : ',i6)
      write(iwr,75) n1, n2, n3, n4, n5
75    format(//,t15,'implicit definition of: ',/,
     &          t39,i4,' internal singles',/    
     &          t39,i4,' external singles',/   
     &          t39,i4,' internal doubles',/     
     &          t39,i4,' mixed    doubles',/    
     &          t39,i4,' external singles',/,   
     &          8x,57(1h*))
      write(iwr,*) 'after removal of redundants'
      write(iwr,*)nconf,nrefco,nsinin,nsinex,ndouin,ndoumi,ndoue1,ndoue2
c...
c...  gather concluding results and write all relevant data to tape15
c...
c...  calculate maxima
      call maxi_crestr(korb,ninner,kvirt,nvirt,kntdet,knstru,nconf,imax,
     &                 maxdet,maxstr)
c...  write
      nrzero = 0
      idoubz(1) = 0
      call tape15(nelec,ndet,nstruc,nconf,kconf,nelec+1,pacdet,kntdet,
     &            knstru,kdetps,idetps,coefff,ncfff,npacd,mult,imax,
     &            maxdet,maxstr,nrzero,idoubz)
c
      write(iwr,85) cpulft(1) - time0
      call end_time_period(TP_VB_STRUC)
85    format(/,t10,'** total cpu time used:      ',f6.2,' cpu-seconds',
     &                                                           ' **')
ccdebug prints           
      call pri(kconf,melec1,idconf,knstru,kntdet,kdetps,idetps,
     &         coefff,pacdet,kdeter,.false.,nconf,nconf)
c
      return
      end
       integer function newtont(i,j)
       implicit REAL (a-h,o-z), integer (i-n)
INCLUDE(../m4/common/iofile)
c...   estimate i over j
       k = i - j
       newtont = 1
       if (j.eq.0.or.j.eq.i) then
          newtont = 1
          return
       end if
       jj = 2
       kk = 2
       do 30 ii=2,i
          newtont = newtont * ii
          if ( mod(newtont,jj).eq.0.and.jj.le.j ) then
             newtont = newtont / jj
             jj = jj + 1
          else if ( mod(newtont,kk).eq.0.and.kk.le.k ) then
             newtont = newtont / kk
             kk = kk + 1              
          end if
30     continue
       if (jj.ne.j+1.or.kk.ne.k+1) then
          write(iwr,11) i,j,newtont 
11        format(' koos is gek volgens newtont ',3i5)
          stop
       end if
       return
       end
      subroutine nobeta (kconf,knwdet,klogic,nnwdet,knwstr,coefco,
     &                   ndetst,nnwstr,melec,mdetst)
c***********************************************************************
c
      implicit REAL  (a-h,o-z) , integer (i-n)
c
c***********************************************************************
c
c     'nobeta' creates both the structure, and the determinant, for the
c     case that there are no single occupied orbitals with alpha-spin.
c     (including the case of a closed shell).  in addition, it assigns
c     values to a few variables, that must be given through to the main
c     program. (these varables are : ndetst, nnwstr and nnwdet).
c
c
      dimension kconf  ( melec+1 )        ,
     &          knwstr ( mdetst  ,   *   ),
     &          coefco ( mdetst  ,   *   ),
     &          knwdet ( melec   ,   *   ),
     &          klogic (*)
c
      common /condep/ nsocc,ndocc,nalpha,nbeta,npair
c
c                           define structure
c
      knwstr(1,1) = 1
      klogic(1)   = 1
      coefco(1,1) = 1
c
      ndetst = 1
      nnwstr = 1
      nnwdet = 1
c
c                          define determinant
c
c                 all single occupieds have alpha-spin
c
      do 10 ialpha=1,nalpha
         knwdet(ialpha,1) = kconf(ialpha)
10    continue
c
c                          double occupieds
c
      do 20 idocc=1,ndocc
         knwdet(nsocc      +idocc,1) = kconf(nsocc+2*idocc)
         knwdet(nsocc+ndocc+idocc,1) = kconf(nsocc+2*idocc)
20    continue
c
      return
      end
      subroutine noproj (nnwstr,knway,kbeta,ndetst,kposar,knwstr,coefco,
     &                   nwval1,nwval2,mpair,mperm,mdetst)
c***********************************************************************
c
      implicit REAL  (a-h,o-z) , integer (i-n)
c
c***********************************************************************
c
c     'noproj' expects the elements '(1,index2,index3,index4)' to
c     '(nbeta,index2,index3,index4)' to be collection of beta-orbitals,
c     belonging to a determinant of structure 'index4'.
c
c     each collection of beta-orbitals is given logic determinant-number
c     based on the orbital-numbers, and a number that is determined by
c     the order of appearance on the array kbeta. the logic determinant-
c     number is used to define structures (including signs) on the array
c     'knwstr'. for each logic determinant-number an appearance-number
c     is put on the array 'kposar'.
c
      dimension knway  ( 0:mpair ),
     &          kbeta  ( mpair   , mperm , 0:mpair , * ),
     &          kposar (   *   ),
     &          knwstr ( mdetst  ,   *   ),
     &          coefco ( mdetst  ,   *   )
c
INCLUDE(common/turtleparam)
c
      common /fac/    kfacul(0:mmsocc)
      common /indep/  time0,nelec,mult
      common /condep/ nsocc,ndocc,nalpha,nbeta,npair
c
      do 10 istruc=1,nnwstr
c
         idet = 0
c
         do 20 inperm=0,npair
            do 30 iway=1,knway(inperm)
c
c            treat determinant 'idet' of structure 'istruc'
c
               idet = idet + 1
c
c              determine logic determinant-number (jlogic)
c
               jlogic = 1
               jcheck = 1
c
               do 110 invest=1,nbeta
c
                  jdiff = kbeta(invest,iway,inperm,istruc) - jcheck
c
                  if (jdiff.eq.0) then
                     jcheck = jcheck + 1
                  else
                     junder = nbeta - invest
                     do 120 ifixed=jcheck,
     &                             kbeta(invest,iway,inperm,istruc)-1
                        jabove = nsocc - ifixed
ckoos                        jlogic = jlogic + kfacul(jabove)/(kfacul(junder)
ckoos     &                                           *kfacul(jabove-junder))
                        jlogic = jlogic + newtont(jabove,junder)
ckoos
120                  continue
                     jcheck = kbeta(invest,iway,inperm,istruc) + 1
                  endif
110            continue
c
c     link the logic determinant-number to an occurance-number of the
c     corresponding beta-orbital-collection on the array 'kbeta'
c
               kposar(jlogic) = (istruc-1)*ndetst + idet
c
30          continue
20       continue
10    continue
c
c                        no spin-projection!!!
c
c     each determinant forms a structure, all coefficients must be set
c     equal to one
c
c     (new) new number of structures : nwval1 (spinprojected : nnwstr)
c     (new) number of determinants per structure : nwval2 (spin-
c                                                  projected : ndetst)
c
ckoos      nwval1 = kfacul(nsocc)/(kfacul(nsocc-nbeta)*kfacul(nbeta))
      nwval1 = newtont(nsocc,nbeta)
ckoos    
      nwval2 = 1
c
      do 200 istruc=1,nwval1
         knwstr(1,istruc) = istruc
         coefco(1,istruc) =      1
200   continue
c
      return
      end
      subroutine normstruc(kdetps,coefff,nstruc)
c...
c...  normalize structures
c...
      implicit REAL (a-h,o-z), integer (i-n)
c
      dimension kdetps(*),
     &          coefff(*)
c
      jcfff  = 1
      do 10 istruc=1,nstruc
c
         cnorm = 0.0d0
         do 20 icfff=jcfff,jcfff-1+kdetps(istruc)
            cnorm = cnorm + coefff(icfff)**2
20       continue
         cnorm = 1.0d0/dsqrt(cnorm)
         do 40 icfff=jcfff,jcfff-1+kdetps(istruc)
            coefff(icfff) = coefff(icfff)*cnorm
40       continue
         jcfff = jcfff + kdetps(istruc)
c
10    continue
c
      return
      end
      subroutine normr(kdetps,idetps,coefff,nostru,nnstru,nncfff)
c...
c...  normalize structures, remove null-functions
c...
      implicit REAL (a-h,o-z), integer (i-n)
c
      dimension kdetps(*),
     &          idetps(*),
     &          coefff(*)
INCLUDE(common/vbcri)
c
      nncfff = 0
      nnstru = 0
      jcfff  = 1
      do 10 istruc=1,nostru
         cnorm = 0.0d0
         do 20 icfff=jcfff,jcfff-1+kdetps(istruc)
            cnorm = cnorm + coefff(icfff)**2
20       continue
         if (cnorm.gt.cridep) then
            if (nnstru+1.lt.istruc) then
               kdetps(nnstru+1) = kdetps(istruc)
               do 30 icfff=1,kdetps(istruc)
                  idetps(nncfff+icfff) = idetps(jcfff-1+icfff)
30             continue
            end if
            cnorm = 1.0d0/dsqrt(cnorm)
            do 40 icfff=1,kdetps(istruc)
               coefff(nncfff+icfff) = coefff(jcfff-1+icfff)*cnorm
40          continue
            nnstru = nnstru + 1
            nncfff = nncfff  + kdetps(istruc)
         end if
         jcfff = jcfff + kdetps(istruc)
10    continue
c
      return
      end
      subroutine occup (kstate,norb)
c***********************************************************************
c
      implicit REAL  (a-h,o-z) , integer (i-n)
c
c***********************************************************************
c
c     'occup' reorders  array 'kstate', by putting
c     double occupied orbitals after single occupieds. the number of
c     double occupied orbitals is stored as row-element 'norb+1'.
c
      dimension kstate ( * )
c
       nsocc = norb
       ndocc =    0
      jorbit =    1
c
10    continue
      jcheck = jorbit + 1
c
20    continue
      if (kstate(jorbit).eq.kstate(jcheck)) then
c
c                put the double twice at the end of the column
c
         jdoub = kstate(jorbit)
c
         do 30 ishift=jorbit,jcheck-2
            kstate(ishift) = kstate(ishift+1)
30       continue
c
         do 40 ishift=jcheck-1,norb-2
            kstate(ishift) = kstate(ishift+2)
40       continue
c
         kstate(norb-1) = jdoub
         kstate(norb  ) = jdoub
c
         nsocc = nsocc - 2
         ndocc = ndocc + 1
      else
         if (jcheck+1.le.nsocc) then
            jcheck = jcheck + 1
            goto 20
         else
            jorbit = jorbit + 1
         endif
      endif
c
      if (jorbit.le.nsocc-1) goto 10
c
c               store the number of double occupied orbitals
c
      kstate(norb+1) = ndocc
c
      return
      end
      subroutine open1(krow,jopen,jdestr)
c...
c...  open up position jopen on cost of the contents of krow(jdestr)
c...
      implicit REAL (a-h,o-z), integer (i-n)
c
      dimension krow(*)
c
      if (jopen.le.jdestr) then
c...  perform right shift
         do 10 ito=jdestr,jopen+1,-1
            krow(ito) = krow(ito-1)
10       continue
      else
c...  perform left shift
         do 20 ito=jdestr,jopen-1, 1
            krow(ito) = krow(ito+1)
20       continue
      end if
c
      return
      end
      subroutine ordernumbs(knumbs,nnumbs)
c...
c...  rearrange the elements of knumbs so that knumbs(i+1) > knumbs(i)
c...  holds for all i. all elements are assumed to be unique.
c...
      implicit REAL (a-h,o-z), integer(i-n)
c
      dimension knumbs(*)
c
      do 10 icycle=1,nnumbs-1
         jcheck = icycle
20       if (knumbs(jcheck).gt.knumbs(jcheck+1)) then
c
                      jaside = knumbs(jcheck  )
            knumbs(jcheck  ) = knumbs(jcheck+1)
            knumbs(jcheck+1) = jaside
            if (jcheck.gt.1) then
               jcheck = jcheck - 1
               go to 20
            end if
c
         end if
10    continue
c
      return
      end
      subroutine orderc(kconf,melec1,nelec,nconf)
c...
c...  order configurations so that in the singles-block 
c...  kconf(iorb+1,iconf) > kconf(iorb,iconf) holds and in the doubles-
c...  block kconf(iorb+2,iconf) > kconf(iorb,iconf) holds.
c...
c...  assumed configuration-construction:
c...       /        singles block        /        doubles block        /
c...               unique numbers                pairs of numbers
c...  size:           nsocc                         2*ndocc  
c...
      implicit REAL (a-h,o-z) , integer(i-n)
c
      dimension kconf(melec1,*)
c
      do 10 iconf=1,nconf
          ndocc = kconf(melec1,iconf)
          nsocc = nelec - 2*ndocc
c...  order singles 
         jfirst = 1
          jlast = nsocc
         do 20 icycle=jfirst,jlast-1
            jcheck = icycle
30          if (kconf(jcheck,iconf).gt.kconf(jcheck+1,iconf)) then
c
                              jaside = kconf(jcheck  ,iconf)
               kconf(jcheck  ,iconf) = kconf(jcheck+1,iconf)
               kconf(jcheck+1,iconf) = jaside
               if (jcheck.gt.jfirst) then
                  jcheck = jcheck - 1
                  go to 30
               endif
c
            endif
20       continue
c...  order doubles
         jfirst = nsocc+1
          jlast = nelec-1
         do 40 icycle=jfirst,jlast-2,2
            jcheck = icycle
50          if (kconf(jcheck,iconf).gt.kconf(jcheck+2,iconf)) then
c
                              jaside = kconf(jcheck  ,iconf)
               kconf(jcheck  ,iconf) = kconf(jcheck+2,iconf)
               kconf(jcheck+1,iconf) = kconf(jcheck+2,iconf)
               kconf(jcheck+2,iconf) = jaside
               kconf(jcheck+3,iconf) = jaside
               if (jcheck.gt.jfirst) then
                  jcheck = jcheck - 2
                  go to 50
               endif
c
            endif
40       continue
c...
10    continue
c
      return
      end
      subroutine orderp(knumbs,nnumbs,nperm)
c...
c...  rearrange the elements of knumbs so that knumbs(i+1) > knumbs(i)
c...  holds for all i. all elements are assumed to be unique.
c...  store the number of permutations, needed to establish the ordering
c...
      implicit REAL (a-h,o-z), integer(i-n)
c
      dimension knumbs(*)
c
      nperm = 0
      do 10 icycle=1,nnumbs-1
         jcheck = icycle
20       if (knumbs(jcheck).gt.knumbs(jcheck+1)) then
c
                      jaside = knumbs(jcheck  )
            knumbs(jcheck  ) = knumbs(jcheck+1)
            knumbs(jcheck+1) = jaside
            nperm = nperm + 1
            if (jcheck.gt.1) then
               jcheck = jcheck - 1
               go to 20
            end if
         end if
c
10    continue
c
      return
      end
      subroutine orders(knumbs,nrow,ncol,jfirst,jlast)
c...
c...  within each column, rearrange the elements of knumbs so that
c...  knumbs(irow+1,icol) > knumbs(irow,icol) holds for 
c...  jfirst <= irow >= jlast. all elements of an ordered group are
c...  assumed to be unique.
c...
      implicit REAL (a-h,o-z), integer(i-n)
c
      dimension knumbs(nrow,*)
c
      do 10 icol=1,ncol
c...
         do 20 icycle=jfirst,jlast-1
            jcheck = icycle
30          if (knumbs(jcheck,icol).gt.knumbs(jcheck+1,icol)) then
c
                              jaside = knumbs(jcheck  ,icol)
               knumbs(jcheck  ,icol) = knumbs(jcheck+1,icol)
               knumbs(jcheck+1,icol) = jaside
               if (jcheck.gt.jfirst) then
                  jcheck = jcheck - 1
                  go to 30
               end if
            end if
c
20       continue
c...
10    continue
c
      return
      end
c      subroutine wrir2d (array,ndim,object,lencol,jfirst,ncol)
cc***********************************************************************
cc
c      implicit REAL  (a-h,o-z) , integer (i-n)
cc
cc***********************************************************************
cc
cc     'wrir2d' writes the two-dimensional real-array 'array' to tape6.
cc     the columns of this array are seen as units (with length 'lencol')
cc     and they are written on one line, preceeded by the object-name
cc     (given through by the calling program by means of character dummy
cc     argument 'object') and the column-number. the columns 'jfirst' to
cc     'ncol' are written.
cc
c      dimension array (ndim, * )
cc
c      character *(*) object
cc
c      write (6) ' '
c      do 10 icol=jfirst,ncol
c         write(iwr,5) object, icol+1-jfirst,
c     &                               (array(irow,icol), irow=1,lencol)
c5       format (' ',a,i4,' : ',20f6.2)
c10    continue
cc
c      return
c      end
      subroutine pair (kterm,nnwstr,kpair,msocc)
c***********************************************************************
c
      implicit REAL  (a-h,o-z) , integer (i-n)
c
c***********************************************************************
c
c     'pair' deals with the first 'nnwstr' columns of the array 'kterm'.
c     these are seen as leading terms. each element of kterm is eighter
c     1 or 0, standing for alpha- and beta-spin respectively. every beta
c     is paired to the closest (unlinked) lower alpha, starting with the
c     lowest beta. in this way 'npair' pairs are formed, npair stemming
c     from common block 'condep'
c
c     for every pair the positions of the linked alpha's and beta's on
c     the array 'kterm' is put on the array 'kpair'. (first alpha's than
c     beta's)
c
      dimension kterm ( msocc , * ),
     &          kpair ( msocc , * )
c
      common /condep/ nsocc,ndocc,nalpha,nbeta,npair
c
c                      treat all leading terms
c
      do 10 iterm=1,nnwstr
c
         jbeta = 2
c
c                link npair beta's to (npair) alpha's
c
         do 20 ipair=1,npair
c
c                find lowest (unlinked) beta-element
c
30          continue
            if (kterm(jbeta,iterm).ne.0) then
c
c                 this is no beta, try next element
c
               jbeta = jbeta + 1
               goto 30
            else
c
c      find nearest unlinked alpha-element (lower than just found beta)
c
               jalpha = jbeta - 1
c
40             continue
               if (kterm(jalpha,iterm).ne.1) then
c
c             this is not an unlinked alpha, try former element
c
                  jalpha = jalpha - 1
                  goto 40
               else
c
c     put respectively on the array kpair : the positions of the just
c                      found alpha and beta on the array kterm
c
                  kpair((ipair-1)*2+1,iterm) = jalpha
                  kpair((ipair-1)*2+2,iterm) = jbeta
c
c                       mark the linked alpha
c
                  kterm (jalpha,iterm) = 2
c
               endif
            endif
c
            jbeta = jbeta + 1
c
20       continue
10    continue
c
      return
      end
      subroutine permut (kpair,knway,nknway,kbeta,nnwstr,
     *                   mpair,mperm,msodd,kperm)
c***********************************************************************
c
      implicit REAL  (a-h,o-z) , integer (i-n)
c
c***********************************************************************
c
c     'permut' interpretates the columns of the array 'kterm' as 'npair'
c     pairs, each column standing for one leading term. ('npair' comming
c     from common block 'condep'). for every pair the lowest element is
c     assumed to be an alpha-orbital, and the other to be a beta-.
c
c     'permut' generates all possible combinations of inside-pair perm-
c     utations. after each combination of permutations only the beta or-
c     bitals are stored, on the four-dimensional array 'kbeta'. the
c     fourth index of this array is formed by the number of the leading
c     term ('iterm'), the third by the number of inside-pair permuta-
c     tations ('inperm') and the second by the number of the current
c     combination of inperm permutations ('iway').
c
c      first a permutation-scheme is formed, whereby the positions of
c     the pairs that must be permuted, is put on the array 'kperm'.
c     the number of combinations of inperm permutations is put on the
c     array 'knway'. both arrays are later used in order to actually
c     apply the permutation-scheme on all leading terms.
c
      logical last
c
      dimension kpair ( msodd   , * ),
     &          kperm ( mpair   , mperm , mpair   ),
     &          knway ( 0:mpair   ),
     &          kbeta ( mpair   , mperm , 0:mpair   , * )
c
      common /condep/ nsocc,ndocc,nalpha,nbeta,npair
INCLUDE(common/turtleparam)
      common    /fac/ kfacul(0:mmsocc)
c
c     the number of ways upon which zero permutations can be
c                                                  carried out equals 1
c
      knway(0) = 1
      nknway   = 1
c
      do 10 inperm=1,npair
c
c     calculate the number of ways upon which 'inperm' permutations
c                                                    can be carried out
c
ckoos         knway(inperm) =
ckoos     &           kfacul(npair)/(kfacul(npair-inperm)*kfacul(inperm))
        knway(inperm) = newtont(npair,inperm)
ckoos
c
         nknway = nknway + knway(inperm)
c
c           define the first way to carry out 'inperm' permutations
c
         do 20 iperm=1,inperm
            kperm(iperm,1,inperm) = iperm
20       continue
c
c         each new way to carry out 'inperm' permutations is based on
c                                             the former way to do this
c
         do 30 iway=2,knway(inperm)
c
            if (kperm(inperm,iway-1,inperm).ne.npair) then
c
c     increase the pair-number of the last permutation, other permu-
c                                                   tations as previous
c
               do 40 iperm=1,inperm-1
                  kperm(iperm,iway  ,inperm) =
     &            kperm(iperm,iway-1,inperm)
40             continue
c
               kperm(inperm,iway  ,inperm) =
     &         kperm(inperm,iway-1,inperm) + 1
c
            else
c
               jshift = inperm - 1
c
50             continue
               if (kperm(jshift+1,iway-1,inperm)-1.eq.
     &             kperm(jshift  ,iway-1,inperm))      then
c
c                    try again, one permutation lower
c
                  jshift = jshift - 1
                  goto 50
               else
c
c                permutations 1 to 'jshift'-1 as previous
c
                  do 60 iperm=1,jshift-1
                     kperm(iperm,iway  ,inperm) =
     &               kperm(iperm,iway-1,inperm)
60                continue
c
c             increase the pair-number of permutation 'jshift'
c
                  kperm(jshift,iway  ,inperm) =
     &            kperm(jshift,iway-1,inperm) + 1
c
c     the pair-numbers of all permutations starting from 'jshift'
c                                                  must form a sequence
c
                  do 70 iperm=jshift+1,inperm
                     kperm(iperm  ,iway,inperm) =
     &               kperm(iperm-1,iway,inperm) + 1
70                continue
c
               endif
            endif
30       continue
10    continue
c
c     perform (virtually) all possible inner-pair permutations on all
c     leading terms, and put the beta-orbitals on an array
c
      do 110 iterm=1,nnwstr
c
c     zero-permutation, copy all beta-orbitals of the array 'kpair', let
c     lower orbital-numbers occur before higher orbital-numbers
c
         do 120 icp=1,npair
            kbeta( icp       ,1,0,iterm) =
     &      kpair((icp-1)*2+2,    iterm)
c
            jcheck = icp
c
130         continue
            if (jcheck.gt.1) then
               if (kbeta(jcheck,1,0,iterm).lt.kbeta(jcheck-1,1,0,iterm))
     &                                                              then
c
c     swap the two elements, so that the lowest orbital-number occurs
c     before the highest orbital-number
c
                  jaside = kbeta(jcheck-1,1,0,iterm)
                  kbeta(jcheck-1,1,0,iterm) = kbeta(jcheck,1,0,iterm)
                  kbeta(jcheck,1,0,iterm) = jaside
c
                  jcheck = jcheck - 1
                  goto 130
               endif
            endif
120      continue
c
         do 140 inperm=1,npair
            do 150 iway=1,knway(inperm)
c
c     get beta-orbitals, so that one of the  'knway(inperm)' possible
c     combinations of 'inperm' inner-pair-permutations is satisfied. for
c     every collection of beta-orbitals : let lower orbital-numbers
c     occur before higher orbital-numbers
c
                last = .false.
               jperm =       1
c
               do 160 ipair=1,npair
                  if (last.or.ipair.ne.kperm(jperm,iway,inperm)) then
c
c         copy the beta-orbital of pair 'ipair' of the array 'kpair'
c
                     kbeta( ipair       ,iway,inperm,iterm) =
     &               kpair((ipair-1)*2+2,            iterm)
                  else
c
c         copy the alpha-orbital of pair 'ipair' of the array 'kpair'
c
                     kbeta( ipair       ,iway,inperm,iterm) =
     &               kpair((ipair-1)*2+1,            iterm)
c
                     if (jperm.eq.inperm) then
                        last = .true.
                     else
                        jperm = jperm +1
                     endif
                  endif
c
                  jcheck = ipair
c
170               continue
                  if (jcheck.gt.1) then
                     if (kbeta(jcheck  ,iway,inperm,iterm).lt.
     &                   kbeta(jcheck-1,iway,inperm,iterm)) then
c
c     swap the two elements to make the lowest orbital-number to occur
c                                   on the array first
c
                        jaside = kbeta(jcheck-1,iway,inperm,iterm)
                        kbeta(jcheck-1,iway,inperm,iterm) =
     &                  kbeta(jcheck  ,iway,inperm,iterm)
                        kbeta(jcheck  ,iway,inperm,iterm) = jaside
c
                        jcheck = jcheck - 1
                        goto 170
                     endif
                  endif
160            continue
150         continue
140      continue
110   continue
c
      return
      end
      subroutine prconf(knstru,kntdet,pacdet,kdetps,idetps,coefff,
     *                  nc,kdeter,melec1,kconf,iconf,ilead)
c
      implicit REAL  (a-h,o-z) , integer (i-n)
c
c
c...    print a configuration
c
INCLUDE(common/c8_16vb)
INCLUDE(common/turtleparam)
INCLUDE(../m4/common/iofile)
      dimension pacdet(*),kdetps(knstru),idetps(*),coefff(*)
      dimension kdeter(melec1,kntdet),kconf(melec1+1)
      common /indep/ time0,nelec,mult
      parameter (maxpr=10000)
      dimension isp(maxpr),isignp(maxpr)
      dimension ilead(*)
      character*60 form
      character*1 ab(2)
      character*1 ldt(melec)
      data ab/'a','b'/
c
      if(nelec.gt.maxpr.or.kntdet.gt.maxpr)call vberr('prconf overflow')
      nsingl = nelec - 2*kconf(melec1+1)
      nalpha = (nelec + mult - 1)/2
c
      nc = 0
c
      call unpack(pacdet,n8_16,kdeter,nelec*kntdet)
c
c ...   translate determinants to spin-pattern corresponding to kconf
c
      do 5 i=1,kntdet
         do 1 k=1,nalpha
1        isp(k) = 1
         do 2 k=nalpha+1,nelec
2        isp(k) = 2
         isignp(i) = ialign(kconf,kdeter(1,i),isp,nelec)
         do 3 k=1,nelec
3        kdeter(k,i) = isp(k)
5     continue
c
c...   structures
c
      write(iwr,600) iconf,(kconf(i),i=1,nelec)
600   format(/' ** configuration ',i5,'  **',/,(t4,15i4))
c      write(iwr,601)
c601   format(/,5x,'******** structures ********',/)
      ioff=0
      do 10 istruc=1,knstru
         do ii=1,nsingl 
            ldt(ii)=ab(ilead(ioff+ii))
         enddo
         ioff=ioff+nelec
         write(iwr,602) istruc,kdetps(istruc),(ldt(ii),ii=1,nsingl)
602      format (5x,'structure',i4,'  determinants',i4,' leading term  '
     *      ,(100a1))
c..              5x,'determinant',2x,'coefficient  ......')
         write(iwr,603) (idetps(idet+nc),coefff(idet+nc)*
     &               isignp(idetps(idet+nc)),idet=1,kdetps(istruc))
603      format (1x,(t5,9(i4,f9.5)))
10    nc = nc + kdetps(istruc)
c
c...    determinants
c
      if (nsingl.le.0) return
      write(iwr,604)
604   format (5x,'******** determinants ********')
c
      nn = 75/(7+nsingl)
      write(form,'(a,i2,a,i2,a)')
     *      '(1x,(t5,',nn,'(i4,3h : ,',nsingl,'a1),2x))'
      write(iwr,form) (idet,(ab(kdeter(i,idet)),
     &                 i=1,nsingl),idet=1,kntdet)
c
      return
      end
      subroutine prepli(jfirst,jlast ,idconf,
     &                  knstru,kntdet,kdetps,pacdet,
     &                  keqcf ,neqcf ,ksocc ,nsocc ,
     &                  jstruc,jcfff ,jpacd ,mstruc,mdet  ,
     &                  kdetp ,idetp ,kcoeff,kdeter,kkdete,idetnr,kcoec,
     &                  ksbeta,nsbeta,kcurcf,mscr)
c...
c...  prepare for call to lindep by remred
c...
      implicit REAL (a-h,o-z), integer (i-n)
c
INCLUDE(common/c8_16vb)
INCLUDE(common/turtleparam)
c
      common /indep/ time0,nelec,mult
      common /fac/   kfacul(0:mmsocc)
c
      dimension idconf(4,*),
c
     &          knstru(*),
     &          kntdet(*),
     &          kdetps(*),
     &          pacdet(*),
c
     &          keqcf (*),
     &          ksocc (*),
     &          kcurcf(*)
c...  gather explicit adresses of equal configurations
      neqcf = 0
      do 10 iconf=jfirst,jlast
         keqcf(neqcf+1) = idconf(4,iconf)
         neqcf = neqcf + 1
10    continue
c...  determine start adresses first configuration
      jstruc = 1
      jcfff  = 1
      jpacd  = 1 
      do 40 iconf=1,keqcf(1)-1
         do 50 istruc=jstruc,jstruc-1+knstru(iconf)
            jcfff = jcfff + kdetps(istruc)
50       continue   
         jstruc = jstruc + knstru(iconf)
         jpacd  = jpacd  + (nelec*kntdet(iconf)-1)/(64/n8_16) + 1
40    continue
c...  gather singly occupieds from first determinant
      call unpack(pacdet(jpacd),n8_16,kcurcf,nelec*kntdet(keqcf(1)))
      nsocc = 0
      do 20 iorb1=1,nelec
         do 30 iorb2=1,nelec
            if (iorb2.eq.iorb1) go to 30
            if (kcurcf(iorb1).eq.kcurcf(iorb2)) go to 20
30       continue
         ksocc(nsocc+1) = kcurcf(iorb1)
         nsocc = nsocc + 1
20    continue
      call ordernumbs(ksocc,nsocc)
c...  calculate mstruc(nsocc) and mdet(nsocc) (conform initsy)
      spin = (mult-1.0d0)/2.0d0
      n2s  = (nsocc/2.0d0-spin+.0001d0)
      if (n2s.eq.0) then
         mstruc = 1
      else
         new1   = newtont(nsocc,n2s)
         n2s1   = n2s-1
         new2   = newtont(nsocc,n2s1)
         mstruc = new1 - new2
      end if
      if (mstruc.gt.maxspvb) call vberr('maxspvb exceeded (prepli)')
      nsalph = (nsocc + mult - 1)/2
      mdet   = newtont(nsocc,nsalph)
      nsbeta = nsocc - nsalph
c...  dynamical core allocation for lindep 
c...  (keqcf and ksocc are on kscr as well)
      kkeqcf = 1
      kksocc = kkeqcf + neqcf
      kdetp  = kksocc + nsocc
      idetp  = kdetp  + mstruc
      kdeter = idetp  + mstruc*mdet
      kkdete = kdeter + mdet*nelec
      idetnr = kkdete + mdet*nelec
      ksbeta = idetnr + mdet
      kendin = ksbeta + nsbeta             
c        
      kcoeff = 1
      kcoec  = kcoeff +   mstruc*mdet     
      kendre = kcoec  + 2*mstruc*mdet
c
      if (kendin.gt.mscr) 
     &   call vberr('core fault in prepli before call to lindep')
c
      return
      end
      subroutine pri(kconf,melec1,idconf,knstru,kntdet,kdetps,
     &               idetps,coefff,pacdet,kdeter,conf,nconf,nnconf)
c...
c...  print vital arrays
c...
      implicit REAL (a-h,o-z), integer (i-n)
c
INCLUDE(common/c8_16vb)
      common /indep/ time0,nelec,mult
INCLUDE(../m4/common/iofile)
c
      dimension kconf (melec1,*),
     &          idconf(4,*),
     &          pacdet(*),
     &          kntdet(*),
     &          knstru(*),
     &          kdetps(*),
     &          idetps(*),
     &          coefff(*),
c
     &          kdeter(melec1-1,*)
c
      logical conf
c
      character*15 object
c
      write(iwr,*) 'vital arrays'
      write(iwr,'(/)')
      nstruc = 0
      ncfff  = 0
      npacd  = 0
      do 10 iconf=1,nconf
         if (conf) then
            write(iwr,*) 'kconf'
            write(iwr,*) (kconf(i,iconf), i=1,nelec)
            print*
         end if
         do 20 iiconf=1,nnconf
            if (idconf(4,iiconf).eq.iconf) then
               write(iwr,*) 'idconf ',(idconf(i,iiconf), i=1,4)
            end if
20       continue
         print*
         write(iwr,*) 'knstru ', knstru(iconf)
         print*
         write(iwr,*) 'kntdet ', kntdet(iconf)
         print*
         write(iwr,*) 'kdetps '
         write(iwr,*) (kdetps(i), i=nstruc+1,nstruc+knstru(iconf))
         print*
         write(iwr,*) 'idetps,coefff'
         do 30 istruc=nstruc+1,nstruc+knstru(iconf)
            write(iwr,*) (idetps(i), i=ncfff+1,ncfff+kdetps(istruc))
            write(iwr,*) (coefff(i), i=ncfff+1,ncfff+kdetps(istruc))
            ncfff = ncfff + kdetps(istruc)
30       continue
         nstruc = nstruc + knstru(iconf)
         print*
         call unpack(pacdet(npacd+1),n8_16,kdeter,kntdet(iconf)*nelec) 
         npacd = npacd + (kntdet(iconf)*nelec-1)/(64/n8_16) + 1 
         write(iwr,*) 'kdeter'
         do 40 idet=1,kntdet(iconf)
            write(iwr,*) (kdeter(i,idet), i=1,nelec)
40       continue
         print*
10       continue
c
      return
      end
      subroutine projc(kconf ,melec1,nconf ,jconf,nspath,ispath,maxsp,
     &                 nofail,
     &                 knstru,kntdet,kdetps,idetps,coefff,pacdet,
     &                 nwstr ,nwcfff,nwpac ,nwdet ,
     &                 krejcf,nrejcf,kstruc,coeff ,kdeter,klogic,
     &                 kscr  ,mscr,scr,irumer)
c...
c...  define structures for nconf configurations
c...
      implicit REAL (a-h,o-z), integer (i-n)
c
INCLUDE(common/c8_16vb)
INCLUDE(common/turtleparam)
c
      dimension pacdet (*),
     &          kconf  (melec1,*),
     &          knstru (*),
     &          kntdet (*),
     &          coefff (*),
     &          kdetps (*),
     &          idetps (*)
      dimension krejcf (*),
     &          kstruc (*),
     &          coeff  (*),
     &          kdeter (melec1-1,*),
     &          klogic (*),
     &          nspath (*),
     &          ispath (maxsp,*),
     &          kscr   (*),
     &          irumer (*),
     &          scr    (*)
c
      logical nofail
c
      common /indep/  time0,nelec,mult
INCLUDE(common/logics)
      common /hans/   msocc,mstrco,mdetco,mdetst,mperm,mstruc
c
      nwcfff = 0 
      nwstr  = 0
      nwpac  = 0
      nwdet  = 0     
      nrejcf = 0
      do 10 iconf=1,nconf    
c...  create structures (based on leading terms)
c.ab
c.ab    set highest
c.ab
            msp = 0
            do 501 isp=1,maxsp
               msp = max(ispath(isp,iconf),msp)
501         continue
c.ab
         call spinef(kconf(1,iconf),nconf,kdeter,klogic,
     &               kntdet(iconf-nrejcf),kstruc,
     &               coeff,kndets,knstru(iconf-nrejcf),allout,projec,
     &               kscr,mscr,mdetst,nelec,irumer,msp)
         if (knstru(iconf-nrejcf).eq.0) then
            if (nofail) then                  
               call vberr(' unexpected failure of spin-projection in pro
     &jc')
            else       
               krejcf(nrejcf+1) = iconf         
               nrejcf = nrejcf + 1
               go to 10
            end if
         end if
c...  produce final results for structures
         call schmic(coeff,mdetst,kndets,knstru(iconf-nrejcf),kstruc,
     &               scr,kntdet(iconf-nrejcf),coefff(nwcfff+1),
     &               kdetps(nwstr+1),idetps(nwcfff+1),ncc,
     &               nspath(iconf),ispath(1,iconf),
     &               rumer(jconf-1+iconf),kdeter,nelec,klogic,irumer)
c...  pack determinants
         call pack(pacdet(nwpac+1),n8_16,kdeter,
     &             nelec*kntdet(iconf-nrejcf))
         nwcfff = nwcfff + ncc
         nwstr  = nwstr  + knstru(iconf-nrejcf)
         nwpac  = nwpac  + (kntdet(iconf-nrejcf)*nelec-1)/(64/n8_16) + 1
         nwdet  = nwdet  + kntdet(iconf-nrejcf)
10    continue
c...  remove redundant determinants, use kdeter for scratch (unpacking)
c.ab      call redund(nelec,nconf-nrejcf,knstru,kntdet,kdetps,idetps,nwpac,
c.ab     &            pacdet,kdeter,nrejd)
c.ab      nwdet = nwdet - nrejd
c
      return
      end
      
      subroutine rdista(istr,nn,max)
c...
c...  read a string numbers (including 'to') 
c...       stop reading if - 'end' is encountered
c...                       - a spintag is encountered
c...                       - # read numbers equals max   
c...  ('i to end' is interpreted as i,...,max)
c...

      implicit REAL (a-h,o-z), integer (i-n)
c
      dimension istr(*)
      logical to
      character*4  itest, ito, ispace, iend, sptag(8)
      character*1  figure(10)
c
      data ito, ispace, iend /'to  ','    ','end '/
      data (sptag(i), i=1,8)
     &     /'rume','lt  ','bd  ','yk  ','gen ','conf','keku','dist'/
      data (figure(i), i=1,10)   
     &     /'1','2','3','4','5','6','7','8','9','0'/     
c
      nn = 0
      to = .false.
      do 10 i=1,max
         istr(i) = 0
10    continue
c
100   call inpa4(itest)
      if (itest.eq.ispace) then
         if (nn.lt.max) then
            call input
            go to 100
         else
            return
         end if
      end if   
c
      if (itest.eq.iend) go to 1000
      if (itest.eq.'core') then
         if (nn.ne.0) call caserr('core not first in  conf')
         call inpi(nc)
         do i=1,nc*2
            istr(i) = (i+1)/2
         end do
         nn = nn + nc*2
         go to 100
      end if
      if (to) then
         if (locatc(figure,10,itest(1:1)).gt.0) then 
            go to 200
         else   
           call vberr('  missing 2nd arg in ..to.. string')
         end if
      end if
      if (locatc(sptag,8,itest).gt.0) go to 2000
      if (itest.eq.ito) then
         if (nn.le.0) then
            call vberr('  missing 1st arg in ..to.. string')
         else
            to = .true.
            go to 100  
         end if
      end if
      if (locatc(figure,10,itest(1:1)).le.0) 
     &   call vberr(' unrecognized symbol in number string') 
c 
      call cinput(jrec,jump,-1)
      nn = nn + 1
      call inpi(istr(nn))
      go to 100
c     
200   nb = istr(nn)
      call cinput(jrec,jump,-1)
      call inpi(ne)
      if (ne.lt.nb) call vberr('  2nd arg in ..to.. string exceeds 1st')
      do 210 i=nb+1,ne
         nn = nn + 1
         istr(nn) = i  
210   continue
      if (nn.gt.max) return
      to = .false.
      go to 100
c
1000  if (to) then        
         nb = istr(nn)
         ne = max
         do 1010 ii=nb+1,ne
            nn = nn + 1
            istr(nn) = ii       
1010     continue
      end if
      return
c
2000  call cinput(jrec,jump,-1)
c
      return
      end
      subroutine rdsymt(isymi,ndim,irconf,isymb,irrep)
c***********************************************************************
c
      implicit REAL  (a-h,o-z) , integer (i-n)
c
c***********************************************************************
c
c     reads symmetry info
c**** only d2h and subgroups thereof ****
c
c***    for use in symdet also
c...   isymi  = 0 : normal entry / read group-name + orbital assignments
c...          > 0 : initialise group  isymi only
c...          < 0 : read groupname but no dimensions
c...   ndim  = dimension of orbital basis
c
c
INCLUDE(common/turtleparam)
c92      character*8 symr,rname,ireff,iblank,isymb,ites
      character*4 symr,rname,ireff,iblank,isymb,irrep,ites
c
      dimension istr (8)
c
      common/vbtable/nirr,mult(8,8),isymr(2,8),jgroup,ifrep,iro(mxorbvb)
      common /tablc/ symr(8),rname(27)
INCLUDE(../m4/common/iofile)
c
c92      data ireff,iblank/'ref',' '/
      data ireff,iblank/'ref ','    '/
c
      isymb = iblank
cccc       mbas =      0
c
      if (ndim.le.0) call err(33)
c
c                     read point-group symbol
c
c92      if (isymi.le.0) call inpa (isymb)
      if (isymi.le.0) call inpa4 (isymb)
c
c                     find point-group number
c
      do 10 igroup=1,8
         nirr = isymr(1,igroup)
         if (isymb.eq.symr(igroup).or.isymi.eq.igroup) then
            jgroup = igroup
            goto 20
         endif
10    continue
c
      call err(1505)
c
20    continue
c....
c.... symbol recognised
c....
      call cinput(jrec,jump,0)
      if(jrec.ge.jump.or.isymi.gt.0) goto 777
c
c                        read irrep. symbol
c
c92      call inpa(ites)
c92      if (ites.eq.ireff) goto 777
      call inpa4(ites)
      irrep=ites
      if (ites.eq.ireff) goto 777
c
c                        find irrep. number
c
      irconf=locatz(rname(isymr(2,jgroup)+1),nirr,ites)
      if(irconf.eq.0) call err(770)
777   continue
      if (isymi.ne.0) return
cccc       nintt=0
cccc       nextt=0
c...
c...    assign orbital-symmetries
c...
      nns = 0
      call izero (ndmim,iro,1)
900   call input
c92      call inpa (ites)
      call inpa4 (ites)
      ifrep = locatz(rname(isymr(2,jgroup)+1),nirr,ites)
      if (ifrep.eq.0) then
         call cinput(jrec,jump,-1)
         call inpi (ifrep)
         if (ifrep.le.0.or.ifrep.gt.nirr) call caserr (' wrong i.r.')
      end if
      call rdist1 (istr,nn)
      do 910 i=1,nn
      iro(istr(i)) = ifrep
910   continue
      nns = nns + nn
      if (nns.lt.ndim) go to 900
c...
c...   check iro
c...
      do 920 i=1,ndim
      if (iro(i).eq.0) go to 950
920   continue
c...
      if (nns.eq.ndim) return
c...
950   write(iwr,960) (iro(i),i=1,ndim)
960   format(' *** error in symmetry assignment *** ',/,' .. iro .. ',/,
     1      (3x,20i3))
      call caserr (' error in symmetry ')
c...
      end
      subroutine reconf(kconf,melec1,nelec,kscr,jconf)
c...
c...  'refer' reads a configuration from kscr and stores it on column 
c...  jconf of the array kconf.
c...  the double occupied orbitals are stored on higher rows than the
c...  single occupied orbitals, and are ordened. the number of double 
c...  occupied orbitals is stored on row nelec+1.   
c...
      implicit REAL (a-h,o-z), integer (i-n)
c
      dimension kconf(melec1,*),
     &          kscr (*)
c...                        read first orbital
      kconf(1,jconf) = kscr(1)
      if (kconf(1,jconf).le.0) call vberr
     &   ('there''s something wrong in the configuration-definition')
c
      nsocc = 1
      ndocc = 0
c...                    read all other orbitals
      do 10 iread=2,nelec
         if (kscr(iread).le.0) call vberr 
     &      ('there''s something wrong in the configuration-definition')
c
c...               check wether this orbital hasn't occured before
         do 20 icheck=1,nsocc
            if (kconf(icheck,jconf).eq.kscr(iread)) then
c...  a double has been found
               do 30 ishift=icheck,nsocc-1
                  kconf(ishift,jconf) = kconf(ishift+1,jconf)
30             continue
               do 40 ishift=nelec-2*ndocc-1,nelec-2
                  kconf(ishift,jconf) = kconf(ishift+2,jconf)
40             continue
               kconf(nelec-1,jconf) = kscr(iread)
               kconf(nelec  ,jconf) = kscr(iread)
               nsocc = nsocc - 1
               ndocc = ndocc + 1
               goto 10
            endif
20       continue
c...  it was a single   
         kconf(nsocc+1,jconf) = kscr(iread)
         nsocc = nsocc + 1
10    continue
c...              store the number of double occupied orbitals
      kconf(nelec+1,jconf) = ndocc
c...           check on more than twice occuring orbitals
      do 60 icheck=1,nsocc
         do 70 idocc=1,ndocc
            if (kconf(icheck,jconf).eq.kconf(nsocc+2*idocc,jconf)) then
               call vberr('mortal sin!!! (pauli-principle is violated)')
            end if
70       continue
60    continue
c...  order double occupieds
      do 80 icount=1,ndocc-1
         jcheck = icount
90       if (kconf(nsocc+2*jcheck  ,jconf).gt.
     &       kconf(nsocc+2*jcheck+1,jconf)) then
                                  jaside = kconf(nsocc+2*jcheck  ,jconf)
           kconf(nsocc+2*jcheck-1,jconf) = kconf(nsocc+2*jcheck+1,jconf)
           kconf(nsocc+2*jcheck  ,jconf) = kconf(nsocc+2*jcheck+1,jconf)
           kconf(nsocc+2*jcheck+1,jconf) = jaside
           kconf(nsocc+2*jcheck+2,jconf) = jaside
            if (jcheck.gt.1) then
               jcheck = jcheck - 1
               go to 90
            end if
         end if
80    continue
c
      return
      end
      subroutine redcon(kconf,kcopy,krejco,melec1,nconf,nelec,nrejco,
     &                  nspath)
c....    
c.... control routine for the removal of multiple defined configurations
c....
      implicit REAL (a-h,o-z), integer (i-n)
c
      dimension kconf (melec1,*),
     &          kcopy (melec1,*),
     &          krejco(*),
     &          nspath(*) 
c
      call icopy (nconf*melec1,kconf,1,kcopy,1)
      call orderc(kcopy,melec1,nelec,nconf)
      call chkcon(kcopy,melec1,nelec,nconf,nrejco,krejco,nspath)
      call ordernumbs(krejco,nrejco)
      call skipi (kconf,melec1,nconf,krejco,nrejco,kconf)
c
      return
      end
c.ab
c.ab May 1999 : subroutine removed, not needed any more
c.ab
c.ab      subroutine redund(nelec,nconf,knstru,kntdet,kdetps,idetps,npredu,
c.ab     &                  pacdet,kdeter,ntred)
c.abc***********************************************************************
c.abc
c.ab      implicit REAL  (a-h,o-z) , integer (i-n)
c.abINCLUDE(common/c8_16vb)
c.abc
c.abc.....remove redundant determinants (those that aren't used in idetps)
c.ab      dimension knstru(*),kntdet(*),idetps(*),pacdet(*),kdeter(nelec,*),
c.ab     &          kdetps(*)
c.ab      msave = 1
c.ab      msav2 = 1
c.ab      it  = 1
c.ab      itold = 1
c.abcccc      it2 = 1
c.abcccc      it3 = 1
c.ab      ntred = 0
c.abc.....
c.ab      do 60 i=1,nconf
c.ab         call unpack(pacdet(itold),n8_16,kdeter,nelec*kntdet(i))
c.abc.....   remove redundant determinants
c.ab         nredun = 0
c.ab         isave  = kntdet(i)
c.abc.....   mm is the starting point in kdetps for the i-th configuration
c.ab         do 50 kk=1,kntdet(i)
c.abc.....   loop over dets. that are checked
c.ab            k = kk - nredun
c.ab            ithere = 0
c.ab            l = 1
c.ab            mm = msav2
c.ab            do 20 m=1,knstru(i)
c.abc.....      loop over the structures of the i-th configuration
c.ab               do 10 n=1,kdetps(mm)
c.abc.....         loop over the dets. of the m-th structure
c.ab                  if (idetps(msave+l-1).eq.k) ithere = 1
c.ab10             l = l + 1
c.ab20          mm = mm + 1
c.ab            ntot = l-1
c.ab            if (ithere.eq.0) then
c.abc.....      remove the redundant determinant
c.ab               nredun = nredun + 1
c.ab               call icopy((kntdet(i)-k)*nelec,kdeter(1,k+1),1,
c.ab     &                    kdeter(1,k),1)
c.ab               kntdet(i) = kntdet(i)-1
c.ab               mm = msav2
c.ab               l = 1
c.ab               do 40 m=1,knstru(i)
c.ab                  do 30 n=1,kdetps(mm)
c.ab                     if (idetps(msave+l-1).gt.k) then
c.ab                        idetps(msave+l-1) = idetps(msave+l-1) - 1
c.ab                     end if
c.ab30                l = l + 1
c.ab40             mm = mm + 1
c.ab            end if
c.ab50       continue
c.ab         call pack(pacdet(it),n8_16,kdeter,nelec*kntdet(i))
c.ab         it  = it + ((kntdet(i)*nelec-1)/(64/n8_16))+1
c.ab         msav2 = msav2 + knstru(i)
c.ab         msave = msave + ntot
c.ab         itold = itold + ((isave*nelec-1)/(64/n8_16)) + 1
c.ab         ntred = ntred + nredun         
c.ab60    continue
c.ab      npredu = it - 1
c.ab      return
c.ab      end
      subroutine remred(idconf,
     &                  knstru,kntdet,kdetps,idetps,coefff,pacdet,
     &                  nconf ,nstruc,ncfff ,npacd ,ndet  ,
     &                  krejcf,kscr  ,scr   ,melec1,mscr)
c...
c...  'remred' removes redundant configurations and structures. 
c...
c...  all configurations are identified by three labels: idconf(1,iconf),
c...  idconf(2,iconf) and idconf(3,iconf). sets of equal configurations
c...  form sequences on idconf. the explicit position of the 
c...  configuration identified by row iconf is idconf(4,iconf).
c...
c...  first a group of equal configurations is searched for, then only
c...  the first one is saved (excstr=.false.) or all structures are 
c...  combined and a new configuration is defined by selecting a linear
c...  independent set.
c...
      implicit REAL (a-h,o-z), integer (i-n)
c
INCLUDE(common/c8_16vb)
INCLUDE(common/turtleparam)
c
      common /indep/ time0,nelec,mult
c
INCLUDE(common/logics)
      common /size/   nrefco,nsinin,nsinex,ndouin,ndoumi,ndoue1,ndoue2
c
      dimension idconf(4,*),
c
     &          knstru(*),
     &          kntdet(*),
     &          kdetps(*),
     &          idetps(*),
     &          coefff(*),
     &          pacdet(*),
c
     &          krejcf(*),
     &          kscr  (*),
     &          scr   (*)
c
      dimension ktype(7)
      equivalence (ktype(1),nrefco),(ktype(2),nsinin),(ktype(3),nsinex),
     &            (ktype(4),ndouin),(ktype(5),ndoumi),(ktype(6),ndoue1),
     &            (ktype(7),ndoue2) 
c
c###### b y p a s s i n g   c o m p i l e r  e r r o r (may 1990) ######
      save
c#######################################################################
      nrejcf = 0
      jfirst = 1
10    continue
c...  find configurations equal to jfirst
      jlast  = jfirst                
20    if (jlast+1.le.nconf) then
         if (idconf(3,jlast+1).eq.idconf(3,jfirst)) then
            if (idconf(2,jlast+1).eq.idconf(2,jfirst)) then
               if (idconf(1,jlast+1).eq.idconf(1,jfirst)) then
c...  configuration jlast+1 equals jfirst
                  jlast = jlast + 1
                  go to 20
               end if
            end if
         end if
      end if
      if (jlast.gt.jfirst) then
         if (excstr) then
c...  the excitation is on structure basis. no configuration can thus a 
c...  priori be expected to define a complete set of structures. combine 
c...  the structures of all the configurations in the group and  
c...  construct a linear independent set.
c
c...  dynamical core allocation for prepli
            keqcf  = 1
            ksocc  = keqcf  + (jlast-jfirst+1)
            kcurcf = ksocc  + nelec
cccc        kend   = kcurcf + nelec     
c...  prepare for lindep
            call prepli(jfirst,jlast,idconf,knstru,kntdet,kdetps,pacdet,
     &                  kscr(keqcf),neqcf,kscr(ksocc),nsocc,jstruc,
     &                  jcfff,jpacd,mstruc,mdet,kdetp,idetp,kcoeff,
     &                  kdeter,kkdete,idetnr,kcoec,ksbeta,nsbeta,
     &                  kscr(kcurcf),mscr)  
            if (nsbeta.eq.0) go to 60
c...  save start adresses of first configuration in group
            jjstru = jstruc
            jjcfff = jcfff
            jjpacd = jpacd
c...  select a set of linear independent structures and define a new 
c...  configuration       
            call lindep(kscr(keqcf),neqcf,kscr(ksocc),nsocc,mdet,mstruc,
     &                  jjstru,jjcfff,jjpacd,knstru,kntdet,kdetps,
     &                  idetps,coefff,pacdet,nnstru,nndet,nncfff,
     &                  kscr(kdetp),kscr(idetp),scr(kcoeff),
     &                  kscr(kdeter),melec1,.true.,kscr(kkdete),
     &                  kscr(idetnr),scr(kcoec),kscr(ksbeta))      
c...  put combined configuration at position of first configuration of
c...  equals group 
            jconf  = idconf(4,jfirst)
c...  size original configuration
            nostru = knstru(jconf)
            nocfff = 0
            do 74 istruc=jstruc,jstruc-1+nostru
               nocfff = nocfff + kdetps(istruc)
74          continue   
            nopacd  = (nelec*kntdet(jconf)-1)/(64/n8_16) + 1
c...  space for old configuration might have to be altered before it  
c...  can be replaced by  the new one
            knstru(jconf) = nnstru
            kntdet(jconf) = nndet
            if (nnstru.le.nostru) then
               incr =  1
            else
               incr = -1
            end if
            call copi(kdetps(jstruc+nostru),kdetps(jstruc+nnstru),
     &                nstruc-(jstruc+nostru)+1,incr)
            call icopy(nnstru,kscr(kdetp),1,kdetps(jstruc),1)
            nstruc = nstruc + (nnstru-nostru)
            if (nncfff.le.nocfff) then
               incr =  1
            else
               incr = -1
            end if
            call copi(idetps(jcfff+nocfff),idetps(jcfff+nncfff),
     &                ncfff-(jcfff+nocfff)+1,incr)
            call icopy(nncfff,kscr(idetp),1,idetps(jcfff),1)
            call copr(coefff(jcfff+nocfff),coefff(jcfff+nncfff),
     &                ncfff-(jcfff+nocfff)+1,incr)
            call fmove(scr(kcoeff),coefff(jcfff),nncfff)
            ncfff = ncfff + (nncfff-nocfff)
            nnpacd = (nelec*nndet-1)/(64/n8_16) + 1
            if (nnpacd.le.nopacd) then
               incr =  1
            else
               incr = -1
            end if
            call copr(pacdet(jpacd+nopacd),pacdet(jpacd+nnpacd),
     &                npacd-(jpacd+nopacd)+1,incr)
            call pack(pacdet(jpacd),n8_16,kscr(kdeter),nelec*nndet)
            npacd = npacd + (nnpacd-nopacd)
60          continue
c...  put the other configurations in the group at the list of rejected
            call icopy(neqcf-1,kscr(keqcf+1),1,krejcf(nrejcf+1),1)
            nrejcf = nrejcf + (neqcf-1)
         else
c...  the excitation is on configuration basis, so all configurations in
c...  the group define the same set of structures. save only the first
c...  one.
            do 50 iconf=jfirst+1,jlast
               krejcf(nrejcf+1) = idconf(4,iconf)
               nrejcf = nrejcf + 1
50          continue
         end if
      end if
c
      jfirst = jlast + 1
      if (jfirst.lt.nconf)  go to 10
c...
c...  remove redundant configurations
c...
      call ordernumbs(krejcf,nrejcf)
      jrejcf = 1
c...  reference adresses
      jostru = 1
      jocfff = 1
      jopacd = 1
c...  new adresses
      jnconf = 1
      jnstru = 1
      jncfff = 1
      jnpacd = 1
      ndet   = 0
c...  configuration type  
      jtype        = 1
      njtype       = ktype(jtype)
      ktype(jtype) = 0
      jjtype       = 0
c...  copy configuration definitions and forget about the redundants
c...  meanwhile, update the number of configurations per type
      do 100 iconf=1,nconf       
         if (jnconf.gt.nconf-nrejcf) then
            if (jtype.lt.7) then
               do 110 itype=jtype+1,7
                  ktype(itype) = 0
110            continue
            end if
            go to 150
         end if
         jjtype = jjtype + 1
         if (jjtype.gt.njtype) then
120         jtype        = jtype + 1
            njtype       = ktype(jtype)  
            if (njtype.eq.0) go to 120
            ktype(jtype) = 0
            jjtype       = 1
         end if
c...  size of refered configuration
         nncfff = 0
         do 130 istruc=jostru,jostru-1+knstru(iconf)
            nncfff = nncfff + kdetps(istruc)
130      continue
         nnpacd  = (nelec*kntdet(iconf)-1)/(64/n8_16) + 1
         if (jrejcf.le.nrejcf) then
            if (krejcf(jrejcf).eq.iconf) then    
c...  this one has been rejected
               jrejcf = jrejcf + 1
               go to 140
            end if
         end if
c...  shift configuration definitions
         call copi(kntdet(iconf) ,kntdet(jnconf),1             ,1)
         call copi(knstru(iconf) ,knstru(jnconf),1             ,1)
         call copi(kdetps(jostru),kdetps(jnstru),knstru(jnconf),1) 
         call copi(idetps(jocfff),idetps(jncfff),nncfff        ,1)
         call copr(coefff(jocfff),coefff(jncfff),nncfff        ,1)
         call copr(pacdet(jopacd),pacdet(jnpacd),nnpacd        ,1)
c...  raise new adresses
         jnstru = jnstru + knstru(jnconf)
         jncfff = jncfff + nncfff
         jnpacd = jnpacd + nnpacd
         ndet   = ndet   + kntdet(jnconf)
         jnconf = jnconf + 1
         ktype(jtype) = ktype(jtype) + 1
140      continue
c...  raise reference adresses
         jostru = jostru + knstru(iconf)
         jocfff = jocfff + nncfff
         jopacd = jopacd + nnpacd 
100   continue
c
150   continue
      nconf  = nconf  - nrejcf
      nstruc = jnstru - 1
      ncfff  = jncfff - 1
      npacd  = jnpacd - 1
c
      return
      end
      subroutine revirt(kvirt,kin,nvirt,kconf,melec1,nelec,jconf,ioldm,
     &                                                            nrej)
c
c     revirt ordens the virtuals from kin and puts them on kvirt. the
c     reference configuration on kconf is given a simular ordering.
c     virtuals that are single occupied in the reference configuraton 
c     appear last. the number of 'single occupied virtuals' is stored at
c     kvirt(nvirt+1)
c     if revirt encounters a virtual that is double occupied in the 
c     reference configuration it will stop execution, unless the oldm-
c     option is used. in that case it just skips the particular virtual
c     and prints a warning mesage. 
c
      implicit REAL (a-h,o-z) , integer (i-n)
c
      dimension kvirt(*),
     &          kin(*),
     &          kconf(melec1,*)
INCLUDE(../m4/common/iofile)
c
      ndocc  = kconf(nelec+1,jconf)
      nsocc  = nelec - 2*ndocc
c
      nsinv = 0
      nrej  = 0
      do 10 ivirt=1,nvirt
c..      read a virtual
         jvirt = kin(ivirt)
c...     check wether this virtual occurs on the reference configuration
         do 20 iorb=1,nelec
            if (jvirt.eq.kconf(iorb,jconf)) then
c...           the virtual occurs in the reference configuration
               if (iorb.gt.nsocc) then
c...              it is doubly occupied
                  if (ioldm.eq.0) then
                     call vberr
     & (' a virtual is doubly occupied in the reference configuration')
                  else
                     if (ioldm.eq.1) write(iwr,5)
5                       format(/10x,'warning!!!!!'//,
     & '-- a virtual is doubly occupied in the reference configuration')
c...                 shrink table of virtuals
                     do 30 ishift=nvirt-nrej-nsinv,nvirt-nrej-1
                        kvirt(ishift) = kvirt(ishift+1)
30                   continue
                     nrej = nrej + 1
                     go to 10
                  end if
               end if
c...           put the 'single occupied virtual' after the others 
               do 40 ishift=nvirt-nrej-nsinv,nvirt-nrej-1
                  kvirt(ishift) = kvirt(ishift+1)
40             continue
c       isave doesn't work for some obscure reason
c              call isave(kvirt(nvirt-nrej-nsinv+1),
c                         kvirt(nvirt-nrej-nsinv),nsinv-1)
               kvirt(nvirt-nrej) = jvirt
c...           reorden the reference configuration
               jaside = kconf(iorb,jconf)
               do 50 ishift=iorb,nsocc-1
                  kconf(ishift,jconf) = kconf(ishift+1,jconf)
50             continue
c                call isave(kconf(iorb+1,jconf),kconf(iorb,jconf),
c     &                                              nsocc-1-iorb)
               kconf(nsocc,jconf) = jaside
c
               nsinv = nsinv + 1
               go to 10
c
            end if
20       continue
c...     the virtual does not occur in the reference configuration
         kvirt(ivirt-nrej-nsinv) = jvirt
10    continue
c
c...  store number of 'single occupied virtuals'
      kvirt(nvirt-nrej+1) = nsinv
c
      return
      end
      subroutine rzero (array,nelem)
c***********************************************************************
c
      implicit REAL  (a-h,o-z) , integer (i-n)
c
c***********************************************************************
c
c     'rzero' supplies the first 'nelem' elements of the real-array
c     'array' with the value 0.
c
      dimension array ( * )
c
      do 10 ielem=1,nelem
         array(ielem) = 0.0
10    continue
c
      return
      end
      subroutine schmic(coef,ndim,mdetst,nstruc,kdstr,
     *                  coec,ntd,coeff,ndstrf,kdstrf,kkc,
     *                  nsel,isel,rumer,
     *                  kdeter,nelec,klogic,irumer)
c
      implicit REAL  (a-h,o-z) , integer (i-n)
c
c
c...   normalise (rumer) or schmidt orthognalise spin-coefficents
c...    metric is 1  // formal first dimension of coef ndim
c...    (at start all structures contain mdetst determinants)
c...    kdstr contains indices of determinants
c...   select some spin-paths if requested
c...   return info in ready to write format in coeff,ndstrf,kdstrf
c...   get rid of not needed determinants
c
INCLUDE(../m4/common/iofile)
INCLUDE(common/vbcri)
      logical rumer
      dimension coef(ndim,nstruc),kdstr(ndim,nstruc),coec(ntd,nstruc),
     *          coeff(*),ndstrf(*),kdstrf(*),isel(nsel),
     *          kdeter(nelec,ntd),klogic(*)
      dimension irumer(*)
c
c...   first gather coefficents in intermediate array coec
c                                               
      do 30 i=1,nstruc
         do 10 j=1,ntd
            coec(j,i) = 0.0d0        
10       continue
         do 20 j=1,mdetst
            jtemp = locat1(klogic,nstruc*mdetst,kdstr(j,i))
            coec(jtemp,i) = coef(j,i)       
20       continue
30    continue             
c
c...   now orthogonalise (if bd-functies) and normalise
c
      do 100 ic=1,nstruc
         if (.not.rumer) then      
c...   orthogonalise
            do 60 jc=1,ic-1
               cc = 0.0d0
               do 40 i=1,ntd
40             cc = cc + coec(i,ic)*coec(i,jc)
               do 50 i=1,ntd
50             coec(i,ic) = coec(i,ic) - cc*coec(i,jc)
60          continue
         end if
c...    normalise    
         cnorm = 0.0d0
         do 70 i=1,ntd
70       cnorm = cnorm + coec(i,ic)**2              
         if (cnorm.lt.cridep) call vberr('normalis. failure schmic')
         cnorm = 1.0d0/dsqrt(cnorm)
         do 80 i=1,ntd
80       coec(i,ic) = coec(i,ic)*cnorm
100   continue
c
c...  now do selection if requested
c
      if (nsel.gt.0) then
         do 120 i=1,nsel
            kk = isel(i)
            if (kk.gt.isel(i)) then
               write(iwr,'(a)') 'WARNING: Unexpected ordering in schmic'
               write(iwr,'(a)') '-> wiping branching/rumer diagram info'
               do k=1,nstruc*nelec
                  irumer(k)=0
               end do
            endif
            do 110 k=1,ntd      
               coec(k,nstruc+i) = coec(k,kk)      
110         continue
            do 115 k=1,nelec
               irumer(nelec*(i-1)+k)=irumer(nelec*(kk-1)+k)
115         continue
120      continue  
         do 130 i=1,nsel
         do 130 k=1,ntd
130      coec(k,i) = coec(k,i+nstruc)
         do 137 i=1,nstruc-nsel
         do 137 k=1,nelec
137      irumer(nelec*(nsel+i-1)+k)=0
         nstruc = nsel      
c...   selected structures are now in front again and nstruc adjusted
c...   see if we can kill off a few determinants
c         i = 0
c140      i = i + 1
c         if (i.le.ntd) then
c            do 150 k=1,nstruc
c150         if (dabs(coec(i,k)).gt.cridep) go to 180
c...    found a determinant without a purpose in life  get rid of it
c.ab
c            print *,'Removing determinant',i
c            ntd = ntd - 1
c            i=i-1
c            do 160 j=i,ntd
c            do 160 k=1,nelec
c160         kdeter(k,j) = kdeter(k,j+1)
c            do 170 j=i,ntd
c            do 170 k=1,nstruc
c170         coec(j,k) = coec(j+1,k)
c180         go to 140
c         end if
      end if
c
c...   prepare output
c
      kkc = 0
      do 200 i=1,nstruc
         ndstrf(i) = 0
         do 190 j=1,ntd
            if (dabs(coec(j,i)).gt.cridep) then
               ndstrf(i) = ndstrf(i) + 1
               kkc = kkc + 1
               coeff(kkc) = coec(j,i)
               kdstrf(kkc) = j
            end if
190      continue
200   continue        
c
      return
      end
      subroutine singc(kconf ,melec1,nrefco,nspath,ispath,maxsp ,
     &                 korb  ,norb  ,idconf,
     &                 jactiv,jvirt ,
     &                 knstru,kntdet,kdetps,idetps,coefff,pacdet,
     &                 nconf ,nstruc,ncfff ,npacd ,ndet  ,
     &                 kstruc,coeff ,kdeter,klogic,parity,kscr,
     &                 kcurcf,mscr  ,scr   ,mconf, irumer)
c...
c...  'singc' performs the single excitation jactiv --> jvirt to
c...  reference configurations. the singles are not stored, but spin 
c...  projected immediately after their creation.   
c...
c...  the reference configurations and the (temporarily stored) excited 
c...  configurations are constructed as follows:
c...       /        singles block        /        doubles block        /
c...               unique numbers                pairs of numbers
c...  size:           nsocc                         2*ndocc  
c...
c...  this leads to the following construction of determinants:
c...       /        alpha block        /         beta block        /
c...       /   singles   /   doubles   /   singles   /   doubles   /
c...  size: nalpha-ndocc     ndocc      nbeta-ndocc      ndocc
c...  (both doubles blocks are the same)
c...
c...  within each block, orb(i+1) >= orb(i) holds (for both  
c...  configurations and determinants). this ordening can only be
c...  disturbed in determinants by excitation to external virtuals, for
c...  these are always put at the last positions of their spin blocks.
c...  
      implicit REAL (a-h,o-z), integer (i-n)
c
INCLUDE(common/c8_16vb)
INCLUDE(common/turtleparam)

c
      common /indep/ time0,nelec,mult
c                                       
INCLUDE(common/logics)
c
      dimension kconf (melec1,*)
      dimension nspath(*),
     &          ispath(maxsp,*),
c
     &          korb  (*),
     &          idconf(*),
c
     &          knstru(*),
     &          kntdet(*),
     &          kdetps(*),
     &          idetps(*),
     &          coefff(*),
     &          pacdet(*),
c
     &          kstruc(*),
     &          coeff (*),
     &          parity(*),
     &          kdeter(melec1-1,*),
     &          klogic(*),
     &          kcurcf(*),
     &          kscr  (*),
     &          irumer(*),
     &          scr   (*)
c
      logical nofail
c
      nalpha = (nelec + mult - 1)/2
      nofail = .false.
c
      do 10 irefco=1,nrefco
c
         ndocc = kconf(nelec+1,irefco)
         nsocc = nelec - 2*ndocc
c...  determine current positions of active and virtual
         jposa = ifind(kconf(1,irefco),nelec,jactiv)
         jposv = ifind(kconf(1,irefco),nelec,jvirt )
c...  leave out excitations that lead to nul-functions
         if (jposa.eq.0.or.jposv.gt.nsocc) go to 10
         call icopy(nelec+1,kconf(1,irefco),1,kcurcf,1)
c...  perform single excitation
         call excitc(jactiv,jvirt ,jposa ,jposv ,nsocc ,ndocc ,
     &               kcurcf)
c...  define structures from new configuration
         call projc(kcurcf,melec1,1,irefco,nspath(irefco),
     &              ispath(1,irefco),maxsp,nofail,knstru(nconf+1),
     &              kntdet(nconf+1),kdetps(nstruc+1),idetps(ncfff+1),
     &              coefff(ncfff+1),pacdet(npacd+1),nnstr,nncfff,nnpac,
     &              nndet,krejcf,nrejcf,kstruc,coeff,kdeter,klogic,
     &              kscr,mscr,scr,irumer)
         if (nrejcf.eq.1) go to 10
c
         if (jvirt.eq.254) then
c...  there is an external virtual in the singles part of an alpha or 
c...  beta block. put it last in its spin block. (disturb the orbital-  
c...  ordening within the determinants).
            call unpack(pacdet(npacd+1),n8_16,
     &                  kdeter,kntdet(nconf+1)*nelec)
            do 110 idet=1,kntdet(nconf+1)
               if (kdeter(nalpha-ndocc,idet).eq.(2**n8_16-2)) then
                  jpose = nalpha - ndocc
                  npose = nalpha
               else
                  jpose = nelec  - ndocc
                  npose = nelec
               end if
               call open1(kdeter(1,idet),npose,jpose)
               kdeter(npose,idet) = jvirt
               parity(idet) = (-1.0d0)**iabs(npose-jpose)
110         continue
            call pack(pacdet(npacd+1),n8_16,kdeter,
     &                kntdet(nconf+1)*nelec)
            do 120 icfff=ncfff+1,ncfff+nncfff
               coefff(icfff) = parity(idetps(icfff))*coefff(icfff)
120         continue
         end if
c
c...  raise configuration counters 
         nconf  = nconf  + 1
         if (nconf.eq.mconf) call dimerr('configurations',mconf)
         ncfff  = ncfff  + nncfff
         nstruc = nstruc + nnstr
         npacd  = npacd  + nnpac
         ndet   = ndet   + nndet
c...  store configuration identification
         call storid(kcurcf,nconf,idconf,korb,norb)
c
10    continue
c
      return
      end
      subroutine singex (kstate,ndim,norb,jcolgr,jcolnw,kexcit,nexcit,
     &                                                         ngener)
c***********************************************************************
c
      implicit REAL  (a-h,o-z) , integer (i-n)
c
c***********************************************************************
c
c     'singex' performs single and double excitations from the ground-
c     state, defined by column 'jcolgr' of the array 'kstate', to vir-
c     tuals , defined by the array 'kexcit'. the ground-state orbitals
c     may eighter be single or double occupied. a single occupied ground
c     state orbital can also occur as a virtual. the way in wich the
c     ground-state orbitals and the virtuals are ordened by respectively
c     'refer' and 'revirt' in combination with the performance of this
c     subroutine garanties that all created states will be unique. the
c     states are stored as columns on 'kstate' with, for each state,
c     double occupieds occuring after the single occupieds. the number
c     of double occupied orbitals is stored, for each state, as row-
c     element 'norb + 1'.
c
      dimension kstate (ndim, * ),
     &          kexcit ( * )
c
      jstate =                jcolnw
      nexref =      kexcit(nexcit+1)
       ndocc = kstate(norb+1,jcolgr)
       nsocc =        norb - 2*ndocc
c
c          single excitations performed on single occupied orbitals
c
      do 10 iorb=1,nsocc
         do 20 ivirt=1,nexcit-nexref
c
c     with these virtual orbitals a state is created, in wich the same
c     orbitals as in the ground-state are doubly occupied
c
            kstate(norb+1,jstate) =      ndocc
c
            do 30 icp=1,iorb-1
               kstate(icp,jstate) = kstate(icp,jcolgr)
30          continue
            kstate(iorb,jstate) = kexcit(ivirt)
            do 40 icp=iorb+1,norb
               kstate(icp,jstate) = kstate(icp,jcolgr)
40          continue
c
            jstate = jstate + 1
c
20       continue
c
         do 50 ivirt=nexcit-nexref+1,nexcit
c
c     these virtual orbitals are single occupied in the ground-state,
c     they give rise to states that have one extra double occupied orbi-
c     tal. (relative to the ground-state).
c
            if (iorb-(nsocc-nexref).eq.ivirt-(nexcit-nexref)) goto 60
c                c
            kstate(norb+1,jstate) =  ndocc + 1
c
            jdelta = 0
            do 70 icp=1,nsocc
               if (icp.eq.iorb.or.icp.eq.nsocc+ivirt-nexcit) then
                  jdelta = jdelta + 1
                  goto 80
               endif
               kstate(icp-jdelta,jstate) = kstate(icp,jcolgr)
80             continue
70          continue
            kstate(nsocc-1,jstate) = kstate(nsocc+ivirt-nexcit,jcolgr)
            kstate(nsocc  ,jstate) = kstate(nsocc+ivirt-nexcit,jcolgr)
            do 90 icp=nsocc+1,norb
               kstate(icp,jstate) = kstate(icp,jcolgr)
90          continue
c
            jstate = jstate + 1
c
60          continue
50       continue
10    continue
c
c          single excitations performed on double occupied orbitals
c
      do 100 iorb=nsocc+1,norb,2
         do 110 ivirt=1,nexcit-nexref
c
c     with these virtuals a state is created, in wich the number of
c     double occupied orbitals has decreased by one
c
            kstate(norb+1,jstate) = ndocc - 1
c
            do 120 icp=1,nsocc
               kstate(icp,jstate) = kstate(icp,jcolgr)
120         continue
            kstate(nsocc+1,jstate) = kexcit(ivirt)
            kstate(nsocc+2,jstate) = kstate(iorb,jcolgr)
            jdelta = 2
            do 130 icp=nsocc+1,norb
               if (icp.eq.iorb.or.icp.eq.iorb+1) then
                  jdelta = 0
                  goto 140
               endif
               kstate(icp+jdelta,jstate) = kstate(icp,jcolgr)
140            continue
130         continue
c
            jstate = jstate + 1
c
110      continue
c
         do 150 ivirt=nexcit-nexref+1,nexcit
c
c     with these virtual orbitals states are created, in which one doub-
c     le occupied orbital is replaced by another
c
            kstate(norb+1,jstate) = ndocc
c
            jdelta = 0
            do 160 icp=1,nsocc
               if (icp.eq.nsocc+ivirt-nexcit) then
                  jdelta = 1
                  goto 170
               endif
               kstate(icp-jdelta,jstate) = kstate(icp,jcolgr)
170            continue
160         continue
            kstate(nsocc  ,jstate) = kstate(iorb,jcolgr)
            kstate(nsocc+1,jstate) = kexcit(ivirt)
            kstate(nsocc+2,jstate) = kexcit(ivirt)
            jdelta = 2
            do 180 icp=nsocc+1,norb
               if (icp.eq.iorb.or.icp.eq.iorb+1) then
                  jdelta = 0
                  goto 190
               endif
               kstate(icp+jdelta,jstate) = kstate(icp,jcolgr)
190            continue
180         continue
c
            jstate = jstate + 1
c
150      continue
100   continue
c
      ngener = jstate - jcolnw
c
      return
      end
      subroutine sings(kconf ,melec1,nrefco,
     &                 korb  ,norb  ,idconf,
     &                 jactiv,jvirt ,
     &                 knstru,kntdet,kdetps,idetps,coefff,pacdet,
     &                 nconf ,nstruc,ncfff ,npacd ,ndet  ,
     &                 kdeter,parity,krejd ,keqdet,kscr  ,scr   ,kcurcf,  
     &                 mscr  ,mconf)
c...
c...  'sings' performs the single excitation jactiv --> jvirt to
c...  reference structures, by applying the following excitation 
c...  operator: e(a-->v) = v(a)+a(a) + v(b)+a(b)  
c...  (v=virtual, a=active, (a)=alpha, (b)=beta, +=cross)
c...
c...  assumed reference configuration-construction:
c...       /        singles block        /        doubles block        /
c...               unique numbers                pairs of numbers
c...  size:           nsocc                         2*ndocc  
c...  within each block, orb(i+1) > orb(i) holds for all i
c...
c...  assumed reference determinant-construction:
c...       /        alpha block        /         beta block        /
c...       /   singles   /   doubles   /   singles   /   doubles   /
c...  size:  nalpha-ndocc     ndocc      nbeta-ndocc      ndocc
c...  (both doubles blocks are the same)
c...  within each block, orb(i+1) > orb(i) holds for all i. this 
c...  ordening might be disturbed if external virtuals are involved. 
c...
c...  external virtuals are put at the last positions in their spin 
c...  blocks
c...
      implicit REAL (a-h,o-z), integer (i-n)
c
INCLUDE(common/c8_16vb)
      common /indep/ time0,nelec,mult
c
      logical null
c
      dimension kconf (melec1,*),
c
     &          korb  (*),
     &          idconf(*),
c
     &          knstru(*),
     &          kntdet(*),
     &          kdetps(*),
     &          idetps(*),
     &          coefff(*),
     &          pacdet(*),
c
     &          kdeter(melec1-1,*),
     &          parity(*),
     &          krejd (*),
     &          keqdet(*),
     &          kscr  (*), 
     &          scr   (*),
     &          kcurcf(*)
c
      nalpha = (nelec + mult - 1)/2
      nbeta  = nelec - nalpha
c...  initialize reference adresses 
      jrstru = 1
      jrcfff = 1
      jrpacd  = 1
      do 10 irefco=1,nrefco
         ndocc  = kconf(nelec+1,irefco)
         nsocc  = nelec - 2*ndocc
c...  size of reference
         nrcfff = 0
         do 20 istruc=jrstru,jrstru-1+knstru(irefco)
            nrcfff = nrcfff + kdetps(istruc)
20       continue
c...  determine current positions of active and virtual
         jposa  = ifind(kconf(1,irefco),nelec,jactiv)
         jposv  = ifind(kconf(1,irefco),nelec,jvirt )
c...  leave out excitations that lead to nul-functions
         if (jposa.eq.0.or.jposv.gt.nsocc) go to 10000
         call unpack(pacdet(jrpacd),n8_16,kdeter,nelec*kntdet(irefco))
c...  start adresses excited structures
         jnconf = nconf  + 1
         jnstru = nstruc + 1
         jncfff = ncfff  + 1
c...  go for it
         call excits(jactiv,jvirt,jposa,jposv,nalpha,nbeta,ndocc,nsocc,
     &               irefco,jrstru,jrcfff,nrcfff,jnconf,jnstru,jncfff,
     &               nncfff,knstru,kntdet,kdetps,idetps,coefff,null,
     &               kdeter,melec1,parity,krejd ,keqdet,kscr,scr,mscr)
         if (null) go to 10000
         if (jvirt.eq.254) then
c...  there is an external virtual in the singles part of an  
c...  alpha or beta block. now, disturb the ordening of the orbitals by 
c...  putting the externals last in their spin blocks. 
            do 1010 idet=1,kntdet(jnconf)
               if (kdeter(nalpha-ndocc,idet).eq.254) then
                  jpose = nalpha - ndocc
                  npose = nalpha
               else
                  jpose = nelec  - ndocc
                  npose = nelec
               end if
               call open1(kdeter(1,idet),npose,jpose)
               kdeter(npose,idet) = jvirt
               parity(idet) = (-1.0d0)**iabs(npose-jpose)
1010        continue
            do 1020 icfff=ncfff+1,ncfff+nncfff
               coefff(icfff) = parity(idetps(icfff))*coefff(icfff)
1020        continue
         end if
c...  write determinants
         call pack(pacdet(npacd+1),n8_16,kdeter,nelec*kntdet(jnconf))
c...  normalize structures
         call normstruc(kdetps(jnstru),coefff(jncfff),knstru(jnconf))
c...  update size parameters
         nconf  = nconf  + 1
         if (nconf.eq.mconf) call dimerr('configurations',mconf)
         nstruc = nstruc + knstru(jnconf)
         ncfff  = ncfff  + nncfff
         npacd  = npacd  + (nelec*kntdet(jnconf)-1)/(64/n8_16) + 1
         ndet   = ndet   + kntdet(jnconf)
c...  store configuration specifications
         call reconf(kcurcf,melec1,nelec,kdeter(1,1),1)
         call ordernumbs(kcurcf,nsocc)
         call storid(kcurcf,nconf,idconf,korb,norb)
10000    continue
c...  raise reference adresses
         jrstru = jrstru + knstru(irefco)
         jrcfff = jrcfff + nrcfff
         jrpacd = jrpacd + (nelec*kntdet(irefco)-1)/(64/n8_16) + 1
10    continue
c
      return
      end    
      subroutine skipi(ifrom,nrow,ncol,kdel,ndel,ito)
c...
c...  move all columns from ifrom to ito, except those that are 
c...  specified by kdel
c...
      implicit REAL (a-h,o-z), integer (i-n)
c
      integer   ifrom(nrow,*), ito(nrow,*)
      dimension kdel (*)
c
      jdel = 1
      newcol = 1
      do 10 icol=1,ncol
         if (icol.ne.kdel(jdel)) then
            do 20 irow=1,nrow
               ito(irow,newcol) = ifrom(irow,icol)
20          continue
            newcol = newcol +1
         else
            jdel = min(jdel+1,ndel)
         endif
10    continue
c
      return
      end
      subroutine skipl(lfrom,nrow,ncol,kdel,ndel,lto)
c...
c...  move all columns from lfrom to lto, except those that are 
c...  specified by kdel
c...
      implicit REAL (a-h,o-z), integer (i-n)
c
      logical   lfrom(nrow,*), lto(nrow,*)
      dimension kdel (*)
c
      jdel = 1
      newcol = 1
      do 10 icol=1,ncol
         if (icol.ne.kdel(jdel)) then
            do 20 irow=1,nrow
               lto(irow,newcol) = lfrom(irow,icol)
20          continue
            newcol = newcol +1
         else
            jdel = min(jdel+1,ndel)
         endif
10    continue
c
      return
      end
      subroutine skipr(rfrom,nrow,ncol,kdel,ndel,rto)
c...
c...  move all columns from rfrom to rto, except those that are 
c...  specified by kdel
c...
      implicit REAL (a-h,o-z), integer (i-n)
c
      REAL  rfrom(nrow,*), rto(nrow,*)
      dimension kdel (*)
c
      jdel = 1
      newcol = 1
      do 10 icol=1,ncol
         if (icol.ne.kdel(jdel)) then
            do 20 irow=1,nrow
               rto(irow,newcol) = rfrom(irow,icol)
20          continue
            newcol = newcol +1
         else
            jdel = min(jdel+1,ndel)
         endif
10    continue
c
      return
      end
      subroutine spinef (kconf,nconf,knwdet,klogic,nnwdet,knwstr,coefco,
     &                   ndetst,nnwstr,allout,projec,
     &                   icr,mscr,mdettt,mmelec,msp,irumer)
c***********************************************************************
c
      implicit REAL  (a-h,o-z) , integer (i-n)
c
c***********************************************************************
c
c      'spinef' creates the prescription for the building of
c     structures, hereby only using the spin-eigenvalue ,the number of
c     electrons (which are given through by the common block 'indep'),
c     and the number of double occupied electrons (given through as
c     element 'nelec+1'  of array 'kconf')
c
c     this prescription is then applied on the configuration as defined
c     bye the array 'kconf', to build the wanted structures (and dets)
c
c     the first 'nnwdet' columns of the two-dimensional array 'knwdet'
c     form the determinants, used in the structures. the first 'nnwstr'
c     columns of the array 'knwstr' form the structures, defined in
c     terms of 'kdeter'-column-numbers. the sign, that is given to each
c     determinant, is put on a corresponding position on the array
c     'coefco'. each structure consists of 'ndetst' determinants per
c     structure.
c
      common/hans/ msocc,mstrco,mdetco,mdetst,mperm,mstruc
c
      logical allout
      logical projec
c
      character *16 object
c
      dimension kconf  (  * )    ,
     &          knwdet (  mmelec   , * ),
     &          klogic ( * ),
     &          knwstr ( mdettt   , * ),
     &          coefco ( mdettt   , * ),
     &          irumer (*),
     &          icr(mscr)
c
      common /indep/  time0,nelec,mult
      common /condep/ nsocc,ndocc,nalpha,nbeta,npair
c
INCLUDE(../m4/common/iofile)
c
INCLUDE(common/turtleparam)
INCLUDE(common/aivb)
c
      common /aimotmp/ aimo_lterms,adbg9,adbg10,vbcrreadnr
      integer aimo_lterms(maxcon,maxact),anof,adbg9,adbg10,vbcrreadnr
c
      REAL anms(maxact)
c
c                          fill common block
c
       ndocc = kconf(nelec+1)
       if ( aivb_set ) ndocc = aivb_ndocc
       nsocc = nelec - 2*ndocc
       nbeta = (nsocc - (mult-1))/2
       nalpha = nsocc - nbeta
       npair = nbeta
c
c
c     return immediately  if the spin-eigen-value does
c     not fit the number of single occupied orbitals.
c
      if (mult-1.gt.nsocc) then
         nnwdet = 0
         nnwstr = 0
         ndetst = 0
         if (allout) then
            write(iwr,5)
5           format (/,' ','** the number of excess alpha-electrons ex',
     &               'ceeds the number of single occupied orbitals **')
         endif
         return
      endif
c
c...    core partitioning for local arrays
c
c    &          kterm  (  msocc   , mstrco ),
c    &          kpair  (  msocc   , mstrco ),
c    &          knway  ( 0:npair ),
c    &          kbeta  (  npair , mperm , 0:npair   , mstrco ),
c    &          kposar ( mdetco )
c    &          kperm  ( mpair  , mperm , mpair )
c
      kterm  = 1
      kpair  = kterm  + nsocc*mstrco
      knway  = kpair  + nsocc*mstrco
      kbeta  = knway  + npair+1
      kposar = kbeta  + npair*mperm*(npair+1)*mstrco
      kperm  = kposar + mdetco
      kend   = kperm  + npair*mperm*npair - 1
      if (kend.gt.mscr) call vberr(' core fault in spinef')
c
c              treat the special case : only alpha-orbitals
c
      if (nbeta.eq.0) then
         call nobeta (kconf,knwdet,klogic,nnwdet,knwstr,coefco,ndetst,
     &                nnwstr,nelec,mdetst)
c
         if ( aivb_set .eqv. .false. ) then
         call stlead(icr(kterm),icr(kpair),irumer,icr(kend),mscr-kend,
     &               nsocc,nnwstr,nconf,'beta')
         else
           aimo_lterm(1) = 1
         end if
         if (allout) then
            write(iwr,6)
            write(iwr,7)
            write(iwr,8)
            write(iwr,9)
            write(iwr,13)
6           format (/,' ',23x,'**** define structures ****')
7           format (/,' ',4x,'structure   1 :')
8           format (' ',24x,'determinant  coefficient')
9           format (' ',24x,'     1           1.0    ')
13          format (/,' ',22x,'**** define determinants ****')
            object = 'determinant     '
            call writ2d (knwdet,mmelec,object,nelec,1,1)
         endif
         return
      endif
c
c                      perform the spin-projection
c
c                         create leading terms
c
      call crterm (icr(kterm),nsocc,nnwstr,msp)

      if ( aivb_set ) then
        anof = 0
        aivb_nlterm = 0
        do i=1,nnwstr
          ilterm = 0
          do ii=1,aimo_elgrps(nsocc,1)
            anms(ii) = 0.0
          end do
          do j=1,nsocc
            aimo_lterms(i,j) = icr(kterm+j-1+(i-1)*nsocc)
            if(allout) then
              write(iwr,'(a,I3)') ' grp   : ',aimo_elgrps(j,1)
              write(iwr,'(a,I3)') ' spin  : ',aimo_elgrps(j,2)
              write(iwr,'(a,I3)') ' kterm : ',aimo_lterms(i,j)
            end if
            if ( aimo_lterms(i,j) .eq. 1 ) then
              anms(aimo_elgrps(j,1)) = anms(aimo_elgrps(j,1)) + 0.5
            else if ( aimo_lterms(i,j) .eq. 0 ) then
              anms(aimo_elgrps(j,1)) = anms(aimo_elgrps(j,1)) - 0.5
            end if
            if(allout)
     1      print *,'========'
          end do

          if(allout)
     1      write(iwr,'(a,10F4.1)') 
     2      ' nms: ',(anms(ii), ii=1,aimo_elgrps(nsocc,1))

          do ii=1,aimo_elgrps(nsocc,1)
            if ( abs(anms(ii)) .eq. aivb_nms(ii) ) ilterm = ilterm + 1
            if (allout) then
              write(iwr,'(2(a,F4.1))')
     1        ' anms: ',anms(ii),
     2        ' _nms: ',aivb_nms(ii)
              write(iwr,'(a,I3)') ' ilterm: ',ilterm
            end if
          end do
          if ( ilterm .eq. aimo_elgrps(nsocc,1) ) then
            do ii=1,aimo_elgrps(nsocc,1)
              if ( aimo_nnd(ii)*0.5 .ne. DABS(anms(ii)) ) goto 143
            end do
            goto 145
143         continue
            do j=1,nsocc
              if ( aimo_elgrps(j,2) .ne. aimo_lterms(i,j) 
     1        .and. aimo_elgrps(j,2) .ne. 2 ) then
                anof = anof + 1
                if (allout) then
                  write(iwr,'(a,I3,a)') 
     1            ' leading term ',i,' not matching '
                  write(iwr,'(a,I3)') ' not found = ',anof
                end if
                goto 144
              end if
            end do
145         continue
            aivb_nlterm = aivb_nlterm + 1
            if ( aivb_nlterm .eq. aivb_inpspp(vbcrreadnr)) then
              aimo_lterm(aivb_nlterm) = i
              if (allout) write(iwr,'(a,I3)') ' found a leading term: ',
     1        aimo_lterm(aivb_nlterm)
              goto 14
            end if
          end if
144       continue
          if (allout) then
            print *,'------------------'
          end if
        end do
        if ( anof .eq. nnwstr ) then
          call vberr(' no matching leading term found in spinef ')
        end if
14      continue
      end if
c
      if (allout) then
c
c                    write the array kterm to output
c
         write(iwr,15)
15       format (/,' ',21x,'**** create leading terms ****')
         object = 'leading term    '
         call writ2d (icr(kterm),nsocc,object,nsocc,1,nnwstr)
      endif
c
c         interpretate leading terms in terms of paired orbitals
c
      call pair (icr(kterm),nnwstr,icr(kpair),nsocc)
c
      if ( aivb_set .eqv. .true. ) then
c
c                    write the array kpair to output
c
        if ( adbg9 .eq. 1 ) then
          write(iwr,25)
          object = 'leading term    '
          call writ2d (icr(kpair),nsocc,object,nbeta*2,1,nnwstr)
        end if
c
      else
c
c                    write the array kpair to output
c
        if ( allout) then
          write(iwr,25)
25        format (/,' ',26x,'**** assign pairs *****',/,' ',
     &                 '(the pairs are standing two by two)')
          object = 'leading term    '
          call writ2d (icr(kpair),nsocc,object,nbeta*2,1,nnwstr)
        end if
      end if
c
c...  store bond patterns in common block leading
c
      if ( aivb_set .eqv. .false. ) then
      call stlead(icr(kterm),icr(kpair),irumer,icr(kend),mscr-kend,
     &            nsocc,nnwstr,nconf,'norm')
      end if
c
c         make all possible combinations of inner-pair-permutations
c
      call permut (icr(kpair),icr(knway),ndetst,icr(kbeta),nnwstr,
     *             npair,mperm,nsocc,icr(kperm))
c
      if ( aivb_set .eqv. .true. ) then
c
c                    write the array kbeta to output
c
        if ( adbg9 .eq. 1 ) then
          write(iwr,35)
          call writ4d (icr(kbeta),nbeta,icr(knway),npair,nnwstr,
     *                 npair,mperm)
        end if
c
      else
c
c                    write the array kbeta to output
c
        if ( allout) then
           write(iwr,35)
35         format (/,' ',25x,'**** permute pairs ****')
c
           call writ4d (icr(kbeta),nbeta,icr(knway),npair,nnwstr,
     *                  npair,mperm)
        end if
      endif
c
c          define structures in terms of logic determinant-numbers
c             
         call structc(nnwstr,icr(knway),icr(kbeta),ndetst,icr(kposar),
     &                knwstr,coefco,npair,mperm,mdetst)
c
      if ( aivb_set .eqv. .true. ) then
c
c     write structures, defined by determinant-numbers and signs, to
c                                                                output
c
        if ( adbg10 .eq. 1 ) then
         write(iwr,45)
c
         do 11 inwstr=1,nnwstr
c
            write(iwr,55) inwstr
            write(iwr,65)
c
            do 21 inwdet=1,ndetst
             write(iwr,75) knwstr(inwdet,inwstr), coefco(inwdet,inwstr)
21          continue
11       continue
        end if

      else
c
c     write structures, defined by determinant-numbers and signs, to
c                                                                output
c
        if ( allout) then
         write(iwr,45)
45       format (/,' ',23x,'**** define structures ****')
c
         do 10 inwstr=1,nnwstr
c
            write(iwr,55) inwstr
55          format (' ',4x,'structure',i4,' :')
            write(iwr,65)
65          format (' ',24x,'determinant',2x,'coefficient')
c
            do 20 inwdet=1,ndetst
             write(iwr,75) knwstr(inwdet,inwstr), coefco(inwdet,inwstr)
75           format (' ',26x,i8,9x,f5.1)
20          continue
10       continue
        end if
      end if
c
c                         define determinants
c

c      print *,'kposar: ',(icr(kposar+i), i=1,mdetco)
c      print *,'ndetst: ',ndetst
c      print *,'kbeta: ',(icr(kbeta+i), i=1,npair*mperm*(npair+1)*mstrco)
c      print *,'kconf: ',(kconf(i), i=1,nelec1)
c      print *,'nnwdet: ',nnwdet
c      print *,'nelec: ',nelec
c      print *,'npair: ',npair
c      print *,'mperm: ',mperm
c      print *,'nnwstr: ',nnwstr
c      do i=1,mstrco
c        print *,'knwstr: ',(knwstr(j,i), j=1,mdett)
c      end do
c      print *,'mdetst: ',mdetst
c      print *,'klogic: ',(klogic(i), i=1,8)
c      print *,'mdetco before deter: ',mdetco
      call deter (icr(kposar),ndetst,icr(kbeta),kconf,knwdet,nnwdet,
     *            nelec,npair,mperm,nnwstr,knwstr,mdetst,klogic)
c
      if (allout) then
c
c                write the array knwdet to output
c
         write(iwr,85)
85       format (/,' ',22x,'**** define determinants ****')
         object = 'determinant     '
         call writ2d (knwdet,mmelec,object,nelec,1,nnwdet)
      endif
c
      return 
      end
      subroutine storid(kconf,nconf,idconf,korb,norb)
c...
c...  establish identification of configuration nconf and store it on 
c...  idconf so that its ordening remains. the identification consists
c...  of three numbers. these numbers (idconf(1,i) to idconf(3,i))
c...  characterize uniquely configuration idconf(4,i).
c...
c...  idconf(1,i) : number of double occupieds  
c...  idconf(2,i) : logical number singles of configuration i
c...  idconf(3,i) : logical number doubles of configuration i
c...  idconf(4,i) : explicit position of configuration i
c...
c...  ordening of idconf:
c...  idconf(1,i) <= idconf(1,i+1)
c...  for idconf(1,i)=idconf(1,i+1) :
c...      idconf(2,i) <= idconf(2,i+1)
c...      for idconf(1,i)=idconf(1,i+1) and idconf(2,i)=idconf(2,i+1) :
c...          idconf(3,i) <= idconf(3,i+1)
c...
      implicit REAL (a-h,o-z), integer (i-n)
c
      common /indep/ time0,nelec,mult
c
      dimension kconf(*),
     &          korb(*),
     &          idconf(4,*)
c
      if (nconf.le.0) call vberr('nconf le 0 in storid')
c...  establish identification of configuration nconf
      ndocc  = kconf(nelec+1)
      nsocc  = nelec - 2*ndocc
      id1cur = ndocc
      id2cur = lognr(kconf(1),nsocc,1,korb,norb)
      id3cur = lognr(kconf(nsocc+1),ndocc,2,korb,norb)
c...  determine proper column for configuration nconf on idconf (jpos)
c...  start comparance at last configuration
      jpos   = nconf - 1
      if (jpos.eq.0) go to 100
c...  compare number of doubles
1     if (id1cur-idconf(1,jpos)) 11,2,100
11    jpos = jpos - 1
      if (jpos.eq.0) then
         go to 100
      else 
         go to 1
      end if
c...  compare logical number singles    
2     if (id2cur-idconf(2,jpos)) 21,3,100
21    jpos = jpos - 1
      if ((jpos.eq.0).or.(id1cur.ne.idconf(1,jpos))) then
         go to 100
      else
         go to 2
      end if
c...  compare logical number doubles
3     if (id3cur-idconf(3,jpos)) 31,100,100
31    jpos = jpos - 1
      if ((jpos.eq.0).or.(id1cur.ne.idconf(1,jpos)).or.
     &                   (id2cur.ne.idconf(2,jpos))) then
         go to 100
      else
         go to 3
      end if
c...
100   jpos = jpos + 1
c...  store identification of configuration nconf at the right column
      call copi(idconf(1,jpos),idconf(1,jpos+1),(nconf-jpos)*4,-1)
      idconf(1,jpos) = id1cur
      idconf(2,jpos) = id2cur
      idconf(3,jpos) = id3cur
      idconf(4,jpos) = nconf 
c
      return
      end
      subroutine structc(nnwstr,knway,kbeta,ndetst,kposar,knwstr,coefco,
     *                   mpair,mperm,mdetst)
c***********************************************************************
c
      implicit REAL  (a-h,o-z) , integer (i-n)
c
c***********************************************************************
c
c     'struct' expects the elements '(1,index2,index3,index4)' to
c     '(nbeta,index2,index3,index4)' to be collection of beta-orbitals,
c     belonging to a determinant of structure 'index4'.
c
c     each collection of beta-orbitals is given logic determinant-number
c     based on the orbital-numbers, and a number that is determined by
c     the order of appearance on the array kbeta. the logic determinant-
c     number is used to define structures (including signs) on the array
c     'knwstr'. for each logic determinant-number an appearance-number
c     is put on the array 'kposar'.
c
c
      dimension knway  ( 0:mpair),
     &          kbeta  ( mpair   , mperm , 0:mpair   , * ),
     &          kposar (   *   ),
     &          knwstr (mdetst ,   *   ),
     &          coefco (mdetst ,   *   )
c
INCLUDE(common/turtleparam)
      common /fac/    kfacul(0:mmsocc)
      common /indep/  time0,nelec,mult
      common /condep/ nsocc,ndocc,nalpha,nbeta,npair
c
      do 10 istruc=1,nnwstr
c
         idet = 0
c
         do 20 inperm=0,npair
            do 30 iway=1,knway(inperm)
c
c            treat determinant 'idet' of structure 'istruc'
c
               idet = idet + 1
c
c              determine logic determinant-number (jlogic)
c
               jlogic = 1
               jcheck = 1
c
               do 110 invest=1,nbeta
c
                  jdiff = kbeta(invest,iway,inperm,istruc) - jcheck
c
                  if (jdiff.eq.0) then
                     jcheck = jcheck + 1
                  else
                     junder = nbeta - invest
                     do 120 ifixed=jcheck,
     &                             kbeta(invest,iway,inperm,istruc)-1
                        jabove = nsocc - ifixed
ckoos                        jlogic = jlogic + kfacul(jabove)/(kfacul(junder)
ckoos     &                                           *kfacul(jabove-junder))
                     jlogic = jlogic + newtont(jabove,junder)
ckoos
120                  continue
                     jcheck = kbeta(invest,iway,inperm,istruc) + 1
                  endif
110            continue
c
c     link the logic determinant-number to an occurance-number of the
c     corresponding beta-orbital-collection on the array 'beta'
c
               kposar(jlogic) = (istruc-1)*ndetst + idet
c
c          add the logic determinant-number to the current structure
c
               knwstr(idet,istruc) = jlogic
c
c     determine the sign, as it is based on the number of permutations
c
               if (inperm.ne.0) then
                  sign = (-1)**inperm
               else
                  sign = +1
               endif
c
c     corrige the sign, for putting the slater-determinant into a alpha-
c     beta blocked form (deteremine the total amount of alpha's standing
c     to the right of beta's)
c
               jralph = 0
               jrbeta = 0
               do 130 ibeta=nbeta,1,-1
                  jralph = jralph + nelec -
     &                   kbeta(ibeta,iway,inperm,istruc) - jrbeta
                  jrbeta = jrbeta + 1
130            continue
c
               if (jralph.ne.0) then
                  sign = sign*(-1)**jralph
               else
                  sign = sign*(+1)
               endif
c
c                      store the sign on an array
c
               coefco(idet,istruc) = sign
c
30          continue
20       continue
10    continue
c
      return
      end
      block data symc
c***********************************************************************
c***********************************************************************
c
      implicit REAL (a-h,o-z), integer (i-n)
INCLUDE(common/turtleparam)
c
      common/crestn/ struct
      character*44 struct
c
      common/vbtable/nirr,mult(8,8),isymr(2,8),jgroup,ifrep,iro(mxorbvb)
c...      nirr : dimension of group            mult : mult-table
c...     isymr : dimension/offset of groups  jgroup : number of group
c...     ifrep : required representation        iro : irreps of orbitals
      common /tablc/ symr(8),rname(27)
c92      character*8 symr, rname
      character*4 symr, rname
c...      symr : names of groups        rname : names of representations
c-----------------------------------------------------------------
      data mult/ 1,2,3,4,5,6,7,8,
     2           2,1,4,3,6,5,8,7,
     3           3,4,1,2,7,8,5,6,
     4           4,3,2,1,8,7,6,5,
     5           5,6,7,8,1,2,3,4,
     6           6,5,8,7,2,1,4,3,
     7           7,8,5,6,3,4,1,2,
     8           8,7,6,5,4,3,2,1  /
c-----------------------------------------------------------------
c c1
      data symr(1),(isymr(i,1),i=1,2) /'c1  ',1,0/
      data rname(1)         /'a   '/
c cs,ci,c2
      data symr(2),(isymr(i,2),i=1,2) /'cs  ',2,1/
      data (rname(i),i=2,3  )  /'a''  ','a'''' '/
      data symr(3),(isymr(i,3),i=1,2) /'ci  ',2,3/
      data (rname(i),i=4,5  )  /'ag  ','au  '/
      data symr(4),(isymr(i,4),i=1,2) /'c2  ',2,5/
      data (rname(i),i=6,7  )  /'a   ','b   '/
c d2,c2v,c2h
      data symr(5),(isymr(i,5),i=1,2) /'d2  ',4,7/
      data (rname(i),i=8,11 )  /'a   ','b1  ','b2  ','b3  '/
      data symr(6),(isymr(i,6),i=1,2) /'c2v ',4,11/
      data (rname(i),i=12,15)  /'a1  ','a2  ','b1  ','b2  '/
      data symr(7),(isymr(i,7),i=1,2) /'c2h ',4,15/
      data (rname(i),i=16,19)  /'ag  ','bg  ','au  ','bu  '/
c d2h
      data symr(8),(isymr(i,8),i=1,2) /'d2h ',8,19/
      data (rname(i),i=20,27)  /'ag  ','b1g ','b2g ','b3g ',
     &                         'au  ','b1u ','b2u ','b3u '/
c
      data struct/' usually on ed7'/
c-----------------------------------------------------------------
      end
      subroutine tape15(nelec,ntdet,ntstru,nconf,kconf,melec1,
     &                  pacdet,kntdet,knstru,kdetps,idetps,coefff,
     &                  ncfff,nwpac,mult,imax,maxdet,maxstr,ndoubmanual,
     &                  idoubmanual)
c...
c...  'tape15' writes all relevant data to tape15, that forms a medium 
c...  for the guidance of turtle
c...  if not atmol ed7 is used to stay in tune with gamess
c...
      implicit REAL (a-h,o-z), integer (i-n)
c
INCLUDE(common/turtleparam)
c
INCLUDE(../m4/common/sizes)
INCLUDE(common/aivb)
c
      dimension kconf (melec1,nconf),
     &          pacdet(nwpac),
     &          kntdet(nconf),
     &          knstru(nconf),
     &          kdetps(ntstru),
     &          idetps(ncfff),
     &          coefff(ncfff),
     &          iarr(mxorbvb),
     &          idoub(mxorbvb),
     &          idoubmanual(mxorbvb)
c
INCLUDE(common/logics)
      common /size/   nrefco,nsinin,nsinex,ndouin,ndoumi,ndoue1,ndoue2
INCLUDE(../m4/common/iofile)
c       
c...  check which orbitals are doubly occupied in all configurations
      do i=1,mxorbvb
         iarr(i)=0
      enddo
      if (.not.manual) then
c...  check how many times an orbital is mentioned in the complete set
        do i=1,nconf
          do j=1,nelec
            iarr(kconf(j,i))=iarr(kconf(j,i))+1
          enddo
        enddo
c...  if an orbital is mentioned 2 times the number of configurations 
c...  it is doubly occupied in all configurations (so nconf*2).
        ndoub=0
        do i=1,mxorbvb
           if (iarr(i).eq.nconf*2) then
              ndoub=ndoub+1
              idoub(ndoub)=i
           endif
        enddo
      end if

      if (manual) then
        ndoub = ndoubmanual
        do i=1,ndoub
          idoub(i) = idoubmanual(i)
        end do
      end if
c
c...  if all orbitals are doubly occupied there are no "active" orbitals
c...  in the subroutine makepsi0. To prevent this we lower the number of
c...  doubly occupied orbitals by one.
      if ((ndoub*2).eq.nelec) then
         write(iwr,*) 'There seem to be only doubly occupied orbitals'
         idoub(ndoub) = 0
         ndoub = ndoub - 1
      end if
      nalpha = (nelec + mult - 1)/2
      if (mrsd) manual = .true.
c
_IF(atmol)
     rewind (15)
cccc      write(iwr,*) 'in tape15'
cccc      write(iwr,*)'  nelec,nalpha,nconf,ntdet,ntstru,maxdet,maxstr,
cccc     &            nwpac,ncfff,imax,manual'
cccc      write(iwr,*)  nelec,nalpha,nconf,ntdet,ntstru,maxdet,maxstr,
cccc     &            nwpac,ncfff,imax,manual
c...                      variables
cccc  voltooien!!!!
      write (15)  nelec,nalpha,nconf,ntdet,ntstru,maxdet,maxstr,
     &            nwpac,ncfff,imax,manual,ndoub
c
cccc  voltooien!!!!
cccc      if (mrsd) write(15) nrefco,nsinin,nsinex,ndouin,ndoumi,ndoue1,
cccc     &                        ndoue2
c...                      number of determinants per structure
      write(15)  kdetps
c...                      determinant-labels
      write(15)  idetps
c...                      corresponding coefficients
      write(15)  coefff
c...                      number of determinants per configuration
      write(15) kntdet
c...                      number of structures per configuration
      write(15) knstru
c...                      packed determinants
      write(15) pacdet
c...                      configurations
      if (.not.(manual.or.mrsd)) write(15) kconf 
c...                      orbitals that are always doubly occupied
      write(15) idoub
c
_ELSE
c
      if ( aivb_set .eqv. .true. ) k7strucold = k7struc
      k7struc = kscra7vb('k7struc',0,'r','p')
      call wr15vb(nelec,nalpha,nconf,ntdet,ntstru,maxdet,maxstr,
     &            nwpac,ncfff,imax,manual,ndoub,'write')
      call wrt3is(kdetps,ntstru,num8)
      call wrt3is(idetps,ncfff,num8)
      call wrt3s(coefff,ncfff,num8)
      call wrt3is(kntdet,nconf,num8)
      call wrt3is(knstru,nconf,num8)
      call wrt3s(pacdet,nwpac,num8)
      if (.not.(manual.or.mrsd)) 
     &   call wrt3is(kconf,melec1*nconf,num8)
      k7struc = kscra7vb('k7struc',0,'r','a')
      k7doub = kscra7vb('k7doub',ndoub,'i','w')
      call wrt3is(idoub,ndoub,num8)
c
_ENDIF
      return
      end
      subroutine wr15vb(nelec,nalpha,nconf,ntdet,ntstru,maxdet,maxstr,
     &                  nwpac,ncfff,imax,manual,ndoub,test)
c
c...  help writing/reading structures file
c
      implicit REAL (a-h,o-z), integer (i-n)
      character*(*) test
      dimension iparam(13)
INCLUDE(../m4/common/iofile)
c
      k7struc = kscra7vb('k7struc',0,'r','n')
      call search(k7struc,num8)
c
      if (test.eq.'write') then
         iparam(1) = nelec
         iparam(2) = nalpha
         iparam(3) = nconf
         iparam(4) = ntdet
         iparam(5) = ntstru
         iparam(6) = maxdet
         iparam(7) = maxstr
         iparam(8) = nwpac
         iparam(9) = ncfff
         iparam(10) = imax
         iparam(11) = manual
         iparam(12) = ndoub
         call wrt3is(iparam,12,num8)
      else if (test.eq.'read') then
         call readis(iparam,12,num8)
         nelec   = iparam(1)
         nalpha  = iparam(2)
         nconf   = iparam(3)
         ntdet   = iparam(4)
         ntstru  = iparam(5)
         maxdet  = iparam(6)
         maxstr  = iparam(7)
         nwpac   = iparam(8)
         ncfff   = iparam(9)
         imax    = iparam(10)
         manual  = iparam(11)
         ndoub   = iparam(12)
      else
         call caserr(' wrong call to wr15vb ')
      end if
c
      return
      end
c***********************************************************************
      subroutine aimonewsymb(nelec,ntdet,ntstru,nconf,kconf,melec1,
     &                  pacdet,kntdet,knstru,kdetps,idetps,coefff,
     &                  ncfff,nwpac,mult,imax,maxdet,maxstr)
      implicit REAL (a-h,o-z), integer (i-n)
c
      parameter (maxcon = 3000)
c
INCLUDE(../m4/common/sizes)
c
      dimension kconf (melec1,nconf),
     &          pacdet(nwpac),
     &          kntdet(nconf),
     &          knstru(nconf),
     &          kdetps(ntstru),
     &          idetps(ncfff),
     &          coefff(ncfff),
     &          iarr(maxorb),
     &          idoub(maxorb)
c
      logical         allout,limout,projec,mrsd,excstr,manual,allrum,
     &                rumer
      common /logics/ allout,limout,projec,mrsd,excstr,manual,allrum,
     &                rumer(maxcon)
      common /size/   nrefco,nsinin,nsinex,ndouin,ndoumi,ndoue1,ndoue2
INCLUDE(../m4/common/iofile)
INCLUDE(common/scra7vb)
c       
c...  check which orbitals are doubly occupied in all configurations
      do i=1,maxorb
         iarr(i)=0
      enddo
      if (.not.manual) then
c...  check how many times an orbital is mentioned in the complete set
         do i=1,nconf
            do j=1,nelec
               iarr(kconf(j,i))=iarr(kconf(j,i))+1
            enddo
         enddo
      end if
c...  if an orbital is mentioned 2 times the number of configurations 
c...  it is doubly occupied in all configurations (so nconf*2).
      ndoub=0
      do i=1,maxorb
         if (iarr(i).eq.nconf*2) then
            ndoub=ndoub+1
            idoub(ndoub)=i
         endif
      enddo
c
c...  if all orbitals are doubly occupied there are no "active" orbitals
c...  in the subroutine makepsi0. To prevent this we lower the number of
c...  doubly occupied orbitals by one.
      if ((ndoub*2).eq.nelec) then
         write(iwr,*) 'There seem to be only doubly occupied orbitals'
         idoub(ndoub) = 0
         ndoub = ndoub - 1
      end if
      nalpha = (nelec + mult - 1)/2
      if (mrsd) manual = .true.

      k7struc = kbl7vb
      call wr15vb(nelec,nalpha,nconf,ntdet,ntstru,maxdet,maxstr,
     &            nwpac,ncfff,imax,manual,ndoub,'write')
      call wrt3is(kdetps,ntstru,num8)
      call wrt3is(idetps,ncfff,num8)
      call wrt3s(coefff,ncfff,num8)
      call wrt3is(kntdet,nconf,num8)
      call wrt3is(knstru,nconf,num8)
      call wrt3s(pacdet,nwpac,num8)
      if (.not.(manual.or.mrsd)) 
     &   call wrt3is(kconf,melec1*nconf,num8)
      k7doub = iposun(num8)
      call wrt3is(idoub,ndoub,num8)
      kbl7vb = iposun(num8)
c 
      return
      end


      subroutine writ1d (karray,object,jfirst,jlast)
c***********************************************************************
c
      implicit REAL  (a-h,o-z) , integer (i-n)
c
c***********************************************************************
c
c     'writ1d' writes the elements 'jfirst' to 'jlast' of the array
c     'karray' to tape6. for each element, the index is given, preceded
c     by a general object-name. this name is given through by the cal-
c     ling program by means of the character dummy argument 'object'.
c
      character *15 object
c
      dimension karray(*)
INCLUDE(../m4/common/iofile)
c
c              write a one-dimensional array to tape6
c
      write(iwr,3)
3     format (' ')
      do 10 ielem=jfirst,jlast
         write(iwr,5) object,ielem,karray(ielem)
5        format (' ',a15,i8,' : ',i4)
10    continue
c
      return
      end
      subroutine writ2d (karray,ndim,object,lencol,jfirst,ncol)
c***********************************************************************
c
      implicit REAL  (a-h,o-z) , integer (i-n)
c
c***********************************************************************
c
c     'writ2d' writes the two-dimensional array 'karray' to tape6. the
c     columns of this array are seen as units (with length 'lencol'),
c     and they are written on one line, preceeded by the object-name
c     (given through by the calling program by means of character dummy
c     argument 'object') and the column-number. the columns 'jfirst' to
c     'ncol' are written.
c
      dimension karray (ndim, * )
INCLUDE(../m4/common/iofile)
c
      character*15 object
c
      write(iwr,3)
3     format(' ')
      do 10 icol=jfirst,ncol
         write(iwr,5) object, icol,
     &                             (karray(irow,icol), irow=1,lencol+1)
5       format (' ',a15,i8,' : ',(t28,20i4))
10    continue
c
      return
      end
      subroutine writ4d (kbeta,jlast1,knway,ind3l,ind4l,mpair,mperm)
c***********************************************************************
c
      implicit REAL  (a-h,o-z) , integer (i-n)
c
c***********************************************************************
c
c     'writ4d' writes the array 'kbeta' to tape6 (output)
c
c
      dimension kbeta ( mpair   , mperm , 0:mpair   , * ),
     &          knway  ( 0:mpair   )
c
      character *25 text1
      character *25 text2
      character *25 text3
      character *26 text4
INCLUDE(../m4/common/iofile)
c
         text1 = '    permutation-scheme'
         text2 = 'number of permutations'
         text3 = '          leading term'
         text4 = 'positions of beta-orbitals'
c
      do 10 index4=1,ind4l
         do 20 index3=0,ind3l
            do 30 index2=1,knway(index3)
c
               write(iwr,3)
3              format(' ') 
               write(iwr,5) text3,index4
               write(iwr,5) text2,index3
               write(iwr,5) text1,index2
5              format (' ',a,' : ',i4)
c
               write(iwr,15) text4, (kbeta(index1,index2,index3,index4),
     &                                                  index1=1,jlast1)
15             format (' ',1x,a,'    : ',20i4)
c
30          continue
20       continue
10    continue
c
      return
      end
      subroutine gam_to_vb(nelec,mult,nbasis)
c
      implicit REAL  (a-h,o-z) , integer (i-n)
c
INCLUDE(../m4/common/sizes)
INCLUDE(../m4/common/infoa)
c...    infoa     
      nelec = ne
      mult = mul 
      nbasis = num
c
      return 
      end

c
      subroutine allconf_drv(kconf,uconf,nconf,rumer,ispath,nspath,
     &                       nelec,mconf,melec1,maxsp)
c
      implicit none
c
INCLUDE(common/turtleparam)
INCLUDE(common/allconf)
INCLUDE(common/splice)
INCLUDE(../m4/common/iofile)
c
      integer iadd,nposs,idoub,icomb,inotd,nelold,ij
      integer i,j,k,l,ncomb_t,npos,nnocc,nconf,inotall,nnot
      integer ifault,nconfmax,fact,ndoub,nrest,iel,ii,jj,norbold
      integer ucomb(maxcon,melec),uposs(maxcon,melec)
      integer uconf(maxcon,melec)
      integer ispath(maxsp,mconf),nspath(mconf),kconf(melec1,mconf)
      logical rumer(mconf)
      integer nelec,melec1,mconf,maxsp,dummy(1)
      character*8 mode
c
      mode = 'allconf  '
c...  how many doubly occ. orbitals we will have at each 'iel' loop
      idoub = 0

c...  how many doubly occ. orbitals there can be
      ndoub = all_nel/2
      nrest = all_nel - ndoub*2

c...  in case we have more electrons than orbitals, we need different starting values
      if ( all_nel .gt. all_norb ) then
        norbold = all_norb
        nelold = all_nel

        all_nel = all_norb
        all_norb = norbold

        idoub = nelold - norbold
        ndoub = nelold/2
        nrest = nelold - ndoub*2
      end if
c
c...  reset nconf (number of conf.)
      nconf = 0
c
c...  start the loop from the number of electrons
c...  decrease by 1 (one elec. less means we have one extra doubly occ. orb.)
c...  and do that until we approach the max. number of doubly occ. orbitals + 'nrest' (if odd)
c...  first loop naturally (for norb > nel) will be all singly occ. orbitals combinations
c...  next loop, will be one doubly and 'nel-1' singly
c...  and so on, and so on
c
      do iel=all_nel,ndoub+nrest,-1
c...  generate combinations 'iel' from number of orbitals set (combinations)
       call gencomb(all_norb,iel,ucomb,ncomb_t,ifault,maxcon,melec,mode,
     &              dummy)

c...  do this only if we go into doubly occ. orbitals mode
       if ( idoub .gt. 0 ) then
c...  generate combinations 'idoub' from 'iel' set (possibilities)
         call gencomb(iel,idoub,uposs,nposs,ifault,maxcon,melec,mode,
     &                dummy)

c...  for each combination ('ncomb_t' in total) we have 'nposs' possibilities
c...  this gives us 'ncomb*nposs' new configurations for given 'idoub'
         do icomb=1,ncomb_t
           do i=1,nposs
             nconf = nconf + 1
             if ( all_ionic .eq. 1 ) then
               ispath(1,nconf) = 0
               nspath(nconf) = 0
             else
               ispath(1,nconf) = 1
               nspath(nconf) = 1
             end if
             rumer(nconf) = all_rumer
c...  'jj' keeps the track of the index in 'uconf' matrix for given 'nconf'
             jj = 1
c...  'inotall' for 'restr' sub-directive. if goes above 1, we exclude given conf.
             inotall = 0
             do j=1,iel
c...  'inotd' checks if it's doubly occ. or not
               inotd = 0
               do k=1,idoub
c...  if given orb. number meant to be doubly occ.
                 if ( j .eq. uposs(i,k) ) then
                   do l=1,all_nrestr
                     if ( ucomb(icomb,j) .eq. all_restr(l)-ncore_i) then
                       inotall = inotall + 1
                     end if
                   end do
c...  check for 'inotall'
                   if ( inotall .gt. 1 ) then
                     nconf = nconf - 1
                     goto 10
                   end if
c...  set the given 'jj' orbital as a doubly occ.
                   uconf(jj,nconf) = ucomb(icomb,j)
                   jj = jj + 1
                   uconf(jj,nconf) = ucomb(icomb,j)
                   jj = jj + 1
c...  it's not a doubly occ.
                 else
                   inotd = inotd + 1
                 end if
               end do
c...  we didn't encounter any doubly occ., so it's singly occ.
               if ( inotd .eq. idoub ) then
                 uconf(jj,nconf) = ucomb(icomb,j)
                 jj = jj + 1
               end if
             end do
c...  copy the configuration to kconf matrix, with proper ordering (doubles at the end)
             call reconf(kconf,melec1,nelec,uconf(1,nconf),nconf)
10           continue
           end do
         end do

       else

c...  if we have configurations with singly occ. orb. only
         do i=1,ncomb_t
           nconf = nconf + 1
           if ( all_covalent .eq. 1 ) then
             ispath(1,nconf) = 0
             nspath(nconf) = 0
           else 
             ispath(1,nconf) = 1
             nspath(nconf) = 1
           end if
           rumer(nconf) = all_rumer
           do j=1,all_norb
             uconf(j,nconf) = ucomb(i,j)
           end do
c...  copy the configuration to kconf matrix, with proper ordering (doubles at the end)
           call reconf(kconf,melec1,nelec,uconf(1,nconf),nconf)
         end do

       end if

c...  increase the 'idoub' for the next loop
       if ( all_noionic .eq. 1 ) return
       idoub = idoub + 1
      end do

      return
      end

c
      subroutine gencomb(n, r, conf, kount, ifault, maxconf, maxel,
     &                   mode, aimo_c)
c
      implicit none
c
c        Algorithm AS 88  Appl. Statist. (1975) Vol.24, No. 3
c
c        When called once, generates all possible combinations
c        from a group of N items.  Each combination (represented in j as       
c        r ordered integers between 1 and n) is processed within allnr.
c
c        Parameters:-
c       
c        n        integer             input:  The size of the group from which
c                                             the combinations are selected.
c
c        r        integer             input:  The size of each comination.
c
c        j        integer array(r)  workspace: Used by allnr to store
c                                              combinations.
c
c        ifault   integer            output:  Fault indicator, equal to:
c                                             0 if 1 le R le N;
c                                             1 otherwise.
c
      integer maxconf, maxel
      integer r, j(r),ii,conf(maxconf,maxel),ifault,i,kount,nmr
      integer n,ip1,l,aimo_c(r),iel1,iel2,kount_c,idbg
      character*8 mode
c
      idbg = 0
      if ( mode(7:8) .eq. '_d' ) idbg = 1

      ifault = 1
      if (r .lt. 1 .or. r .gt. n) return
      ifault = 0
      kount = 0
      kount_c = 0
      nmr = n - r

c
c        Initialize J(1) to lower limit separately, since lower limit for
c        each index depends on lower limit for previous index
c
      i = 1
      j(1) = 1
c
c        Initialize indices for loops i=1,...,r to lower limits
c
    1 if (i .eq. r) goto 3
      ip1 = i + 1
      do 2 l = ip1, r
    2 j(l) = j(l - 1) + 1
c
c        Update the count (kount) of combinations and process the current
c        combination.  The call to Subroutine job may be replaced by
c        statements to process the current combination.
c
    3 kount = kount + 1
c
c        Copy the working j() array to storage conf() matrix
c
       if ( mode(1:7) .eq. 'allconf' ) then
         do ii=1,r
           conf(kount,ii) = j(ii)
         end do
       else if ( mode(1:6) .eq. 'aimo_c' ) then
         if ( idbg .gt. 0 ) print *,'aimo_c mode'
         do iel1=1,r-1
           do iel2=iel1+1,r
             if ( aimo_c(j(iel1)) .eq. aimo_c(j(iel2)) ) then
               if ( idbg .gt. 0 ) then
                 print *,'j(iel1): ',j(iel1)
                 print *,'j(iel2): ',j(iel2)
                 print *,'grp iel1: ',aimo_c(j(iel1))
                 print *,'grp iel2: ',aimo_c(j(iel2))
                 print *,'skipping: ',(j(ii), ii=1,r)
               end if
               goto 12
             end if
           end do
         end do
         kount_c = kount_c + 1
         do ii=1,r
           conf(kount_c,ii) = j(ii)
         end do
         if ( idbg .gt. 0 ) 
     1     print *,'storing: ',kount_c,'  ',(conf(kount_c,ii), ii=1,r)
12       continue
         if ( idbg .gt. 0 ) 
     1     print *,'---------------------------------------------------'
       else if ( mode(1:6) .eq. 'aimo_d' ) then
         if ( idbg .gt. 0 ) 
     1     print *,'aimo_d mode'
         do ii=1,r
           conf(kount,ii) = j(ii)
         end do
       else
         if ( idbg .gt. 0 ) 
     1     print *,'else mode'
       end if
c
c        Increment the first possible index (of loop i) among indices of
c        loops R, R-1,...,1
c
      i = r
    4 if (j(i) .lt. nmr + i) goto 5
      i = i - 1
c
c        Return after all indices have achieved their upper limits
c
      if (i .le. 0) then
       if ( mode(1:6) .eq. 'aimo_c' ) kount = kount_c
       return
      end if
       goto 4
    5 j(i) = j(i) + 1
       goto 1

      end

      recursive function fact(x) result(answer)
      implicit none
c
      integer, intent(in) :: x
      integer :: answer

      if ( x>= 1 ) then
        answer = x*fact(x-1)
      else
        answer = 1
      end if

      end function fact

c
      subroutine remdts(kdeter,nelec,kndet,idetps,idetst,knstru,
     *                  iperm,iiperm)
c***********************************************************************
c
      implicit REAL  (a-h,o-z) , integer (i-n)
c
c***********************************************************************
c
c     'remdts' removes determinants that are generated, but are not needed.
c     kdeter - array containing the determinants
c     nelec - number of electrons
c     kndet - total number of determinants for this configuration
c     idetps - indices of determinants for each structure
c     idetst - number of determinants per structure
c     knstru - number of structures per configuration
c
      dimension kdeter(nelec,kndet)
      dimension idetps(*),idetst(knstru)
      dimension iperm(kndet),iiperm(kndet)
c
c...  Initialize array for checking presence of determinants in structures   
c
      do 10 i=1,kndet
        iperm(i)=0
   10 continue
c
c...  Create permutation array and count number of different determinants
      kk = 0
      do 30 i=1,knstru
        do 20 j=1,idetst(i)
          kk = kk + 1
      iperm(idetps(kk))=1
   20   continue
   30 continue
      k=0
      do i=1,kndet
        if (iperm(i).eq.1) then
      k=k+1
      iperm(i)=k
      endif
      enddo
c... Create inverse permutation array
      do 40 i=1,kndet
        if (iperm(i).ne.0) then
          iiperm(iperm(i))=i
      endif
   40 continue
c... Decrease number of removals from total number of determinants and
c... perform the permutations in the determinants matrix (kdeter)
      kndet=k
      do 60 i=1,kndet
        if(iiperm(i).ne.i) then
      do 50 j=1,nelec
        kdeter(j,i)=kdeter(j,iiperm(i))
   50     continue
        endif
   60 continue
c... Correct the indices in the idetps array as well
      kk = 0
      do 80 i=1,knstru
        do 70 j=1,idetst(i)
          kk = kk + 1
      idetps(kk)=iperm(idetps(kk))
   70   continue
   80 continue
      end
      subroutine stlead(ilead,ipair,irumer,lconf,mscr,nsoc,nstruc,
     &                  nconf,mode)
c...
c...  Structure leading term and bonding arrangement generating subroutine
c...
c...  Upon leaving this routine 'irumer' will contain the bonding arrangements
c...  and leading term for each structure. This information will be used to
c...  print out leading terms in subroutine prconf and to print out the
c...  bonding arrangements at the end of the turtle program.
c...
      implicit none
      common /indep/  time0,nelec,mult
      REAL time0
      integer nelec,mult,kscra7vb,k7conf
INCLUDE(common/turtleparam)
      integer taken(melec)
c
      integer nsoc,ilead(nsoc,*),ipair(nsoc,*),lconf(nelec+1,*),
     &        nstruc,ioff,j,i,irumer(*),iconf,mscr,nconf,npair,ict
      character*4 mode
c
      save iconf
c
INCLUDE(../m4/common/iofile)
      if ((nelec+1)*nstruc.gt.mscr) call corfait ((nelec+1)*nstruc,mscr,
     & 'in stlead')
      k7conf = kscra7vb('k7conf',nconf*(nelec+1),'i','r')
      call readis(lconf,nconf*(nelec+1),num8)
      ioff=0
10    ioff=ioff+1
      if (irumer(ioff).ne.0) goto 10
      ioff=ioff-1
      if (ioff.eq.0) iconf=0
      iconf=iconf+1
      npair=nelec-lconf(nelec+1,iconf)*2-mult+1
      if (mode.eq.'norm') then
        do i=1,nstruc
          ict=0
          do j=1,nelec
            taken(j)=0
          enddo
          if (lconf(nelec+1,iconf).gt.0) then
            do j=nelec-(lconf(nelec+1,iconf)*2)+1,nelec
              ict=ict+1
              irumer(ioff+ict)=lconf(j,iconf)
              taken(j)=1
            enddo
          end if
          do j=1,npair
            ict=ict+1
            irumer(ioff+ict)=lconf(ipair(j,i),iconf)
            taken(ipair(j,i))=1
          enddo
          do j=1,nelec
            if (taken(j).eq.0) then
              ict=ict+1
              irumer(ioff+ict)=lconf(j,iconf)
            endif
          enddo
          do j=1,nsoc
            if (ilead(j,i).ne.0) then
               irumer(ioff+j)=-1*irumer(ioff+j)
            endif
          enddo
          ioff=ioff+nelec
        enddo
      else if (mode.eq.'beta') then
        ict=0
        if (lconf(nelec+1,iconf).gt.0) then
          do j=nelec-(lconf(nelec+1,iconf)*2)+1,nelec
            ict=ict+1
            irumer(ioff+ict)=-1*lconf(j,iconf)
          enddo
        endif
        do j=1,nelec-lconf(nelec+1,iconf)*2
          ict=ict+1
          irumer(ioff+ict)=lconf(j,iconf)
        enddo
        ioff=ioff+nelec
      else
        call vberr('incorrect mode in stlead')
      endif
      return
      end
      subroutine leadsup(iin,iout,num)
c...
c...  'Leading term support routine'
c...  The array 'iin' contains the bonding arrangements per structure. This array is
c...  huge. For every structure there is also a leading term. The sign of the
c...  the value in 'iin' specifies the spin:
c...
c...  negative --> 1 --> later used as 'a' in subroutine prconf
c...  positive --> 2 --> later used as 'b' in subroutine prconf
c...
c...  The sign is no longer needed, so 'iin' will contain no more negative numbers on
c...  exit.
c...
c...  NB The bonding arrangement is not affected nor used in this routine!
c...
      implicit none
c
      integer num,iin(num),iout(num),i
c
      do i=1,num
         if (iin(i).lt.0) then
            iout(i)=1
         else
            iout(i)=2
         end if
      enddo
      return
      end

************************************************************************
      subroutine kekule(v,kconf,kscr,mconf,melec1,ncore_i,nelec,ispath,
     &                  nspath,numkek,maxsp,nconf,fdist)
************************************************************************
c...
c...  get kekule structures from subroutine makekekule and generate
c...  configurations according to these kekule structures. numkek is 
c...  the no. of orbitals involved in kekule structures/configurations.
c...

      implicit REAL (a-h,o-z), integer (i-n)

INCLUDE(../m4/common/iofile)
INCLUDE(common/turtleparam)
INCLUDE(common/tractlt)
INCLUDE(common/logics)

      dimension kconf(melec1,mconf),
     &          ispath(maxsp,mconf),
     &          nspath(mconf),
     &          kscr(*),
     &          v(*)

      memavl = mconf - melec1
      memreq = numkek*3+numkek/2+nbasis

      if (memreq.gt.memavl) then
         write(iwr,*) ' memory available', memavl, 
     +                ' memory requied', memreq
         call caserr(' memory ran out in subroutine kekule ')
      end if

      kneib  = melec1 + 1
      ineib  = kneib  + numkek + numkek/2
      kmo    = ineib  + numkek
      kao    = kmo    + numkek

      call neiblist(v,kscr(kneib),kscr(ineib),kscr(kmo),kscr(kao),
     +              ncore_i,nelec,numkek,kscr,fdist)
*
*     kao is now free
*
      klist   = kmo    + numkek
      kkekule = klist  + numkek
      memlft  = memavl - numkek*4 - numkek/2 
*
      nneib = 0
      nbond = 0
      katom = 1
      nkekule = 0
*
      call makekekule(kscr(kneib),kscr(ineib),kscr(klist),kscr(kkekule),
     +                nkekule,numkek,nneib,katom,nbond,memlft)
*
      if (nkekule.gt.0) then

         write(iwr,10)
10       format(/,5x,' list of all kekule structures (numbers refer to',
     +               ' atom number in the final geometry)')

         do nn = 1, nkekule
            write(iwr,20) nn,(kscr(kmo-1+kscr(kkekule-1+(nn-1)*
     +                        numkek+k)),k=1,numkek)
20          format(/,i10,' : ',(t15,20i4))
         end do

         write(iwr,30)
30       format(/,5x,' corresponding molecular orbitals',
     +               ' (including core)')

         do nn = 1, nkekule
            write(iwr,40) nn,(kscr(kmo-1+kscr(kkekule-1+(nn-1)*
     +                             numkek+k))+kscr(nelec-numkek+1) +
     +                             ncore_i-1,k=1,numkek)
40          format(/,i10,' : ',(t15,20i4))
         end do
c...
c...  translate 'nkekule' kekule structures to 'nkekule' configurations
c...
         do jconf = 1, nkekule
            do i = 1, numkek
               kscr(nelec-numkek+i) = kscr(kmo-1+kscr(kkekule-1+
     +                                (jconf-1)*numkek+i)) +
     +                                kscr(nelec-numkek+1)-1
            end do
            rumer(nconf-1+jconf) = .true.
            ispath(1,nconf-1+jconf) = 1
            nspath(nconf-1+jconf) = 1
            call reconf(kconf,melec1,nelec,kscr,nconf-1+jconf)
         end do
      else
         write(iwr,50)
50       format(/,5x,' NO kekule structures is found ')
      end if

      nconf = nconf - 1 + nkekule

      return
      end

************************************************************************
      subroutine neiblist(v,kneib,ineib,kmo,kao,ncore_i,nelec,numkek,
     &                    kscr,cridst)
************************************************************************
c...
c...  This subroutine gets coordinates of the molecule from GAMESS and
c...  makes  a list of non-redundant neighbours for each atom
c...
      implicit REAL (a-h,o-z), integer (i-n)

INCLUDE(../m4/common/sizes)
INCLUDE(../m4/common/infoa)
INCLUDE(../m4/common/iofile)
INCLUDE(../m4/common/nshel)
INCLUDE(common/turtleparam)
INCLUDE(common/tractlt)
INCLUDE(common/logics)

      dimension  x(numkek), y(numkek), z(numkek)

      dimension kneib(numkek+numkek/2),
     &          ineib(numkek),
     &          kao(*),
     &          kmo(*),
     &          kscr(*),
     &          v(*)
c...
c...  check which a.o. belongs to which atom (e.g., 1,2,3,etc.) and
c...  store the atom label in klab(*)
c...
      n = 0
      do ii = 1,nshell
         do i = kmin(ii),kmax(ii)
            n = n+1
            kao(n) = katom(ii)
         end do
      end do
c...
      cmcoef = 0.0d0
      cmcri = 1.0d-4

      do m = 1, numkek
         j = kscr(nelec-numkek+m)+ncore_i
c...
c...  check the a.o. with the largest coefficient belongs to which atom
c...  corresponding m.o. would also belong to that atom
c...
         cmcoef = dabs(v((j-1)*nbasis+1))
         kk = 1
         do i = 2, nbasis
            if (dabs(v((j-1)*nbasis+i)).gt.cmcoef) then
               kk = i
               cmcoef = dabs(v((j-1)*nbasis+i))
            end if
         end do
         if (cmcoef.lt.cmcri) 
     +      call vberr('vectors are not available to make kekule')
c...
c...  c(1,j),c(2,j) and c(3,j) are the x, y, and z coordinates, respectively
c...  of jth atom
c...
         kmo(m) = kao(kk)
         x(m) = c(1,kao(kk))
         y(m) = c(2,kao(kk))
         z(m) = c(3,kao(kk))
      end do

      if ((numkek/2)*2.ne.numkek) 
     &      call vberr('no kekule structure is possible')
c...
c...  make a list of neighbours on the base of distance and store it in
c...  kneib(*) for each atom
c...
      if (cridst.lt.1.4d0.or.cridst.gt.3.5d0) cridst = 3.1d0
      ncount = 0
      do i = 1, numkek - 1
         nneib = 0
         do j = i + 1, numkek
            dist = dsqrt((x(i) - x(j))**2 + (y(i) - y(j))**2
     +                            + (z(i) - z(j))**2)
            if (dist.le.cridst) then
               ncount = ncount + 1
               nneib = nneib + 1
               kneib(ncount) = j
               if (nneib.gt.3) 
     +         call caserr('too many neighbours in neiblist')
            end if
         end do
         ineib(i) = nneib
      end do
      ineib(numkek) = 0

      if (allout) then
         write(iwr,10)
10       format(/,10x,' NEIGHBOURS LIST FOR ATOMS')
         write(iwr,20)
20       format(/,14x,'atom  |  non-redundant neighbours')
         n = 0
         do i = 1, numkek
            if (ineib(i).ne.0) then
               write(iwr,30) i,(kneib(j),j= n+1,n+ineib(i))
30             format(/,12x,i4,4x,250i4)
               n = n + ineib(i)
            else
               write(iwr,40) i
40             format(/,12x,i4,5x,' no neib')
            end if
         end do
      end if

      return
      end

************************************************************************
      recursive subroutine makekekule(kneib,ineib,klist,kekule,nkekule,
     &                                numkek,nneib,katom,nbond,memlft)
************************************************************************
*      
      implicit REAL (a-h,o-z), integer (i-n)
*
      dimension kneib(numkek+numkek/2),ineib(numkek)
*
      dimension klist(nbond+2)
      dimension kekule(memlft)
*
      logical present
*
*     katom and each of its non-redundant neighbours is a (double) bond.
*     we just want to make possible combinations of these bonds with the 
*     only condition that an atom can not be present twice (i.e., an
*     atom can not be assinged a double bond twice) in any combination. 
*     whenever we get a combination that has all atoms in it that will 
*     be a kekule structure. 

10    present = .false.

*     check if katom is already present in the candidate structure.

      do mm = 2, nbond, 2
         if (katom.eq.klist(mm)) then
            present = .true.
            exit
         end if
      end do
      if (.not.present) then
         if (ineib(katom).ne.0) then
            do i = 1, ineib(katom)

*              check if ith neighbour of katom is already present in the
*              candidate structure.

               do mm = 2, nbond, 2
                  if (kneib(nneib+i).eq.klist(mm)) then
                     present = .true.
                     exit
                  end if
               end do
               if (.not.present) then

*                 store this double bond in the candidate structure 

                  klist(nbond+1) = katom
                  klist(nbond+2) = kneib(nneib+i)

*                 make a recursion.

                  call makekekule (kneib,ineib,klist,kekule,nkekule,
     +                             numkek,nneib+ineib(katom),katom+1,
     +                             nbond+2,memlft)

*
*              else ith neighbour of katom already has a double bond in the
*              candidate structure so we just hit the end do for the increment
*              in i or returned to previous recursive step.

               end if
               present = .false.
            end do

*        else katom is neither present in the candidate structure nor it has
*        any non-redundant neighbour in the list. this means it has been left
*        alone in this candidate structure so we retrun to the previous 
*        recursion.

         end if
      else if (katom.lt.numkek) then

*        katom is already present in the candidate structure. select the next
*        atom from the list as katom.

         nneib = nneib+ineib(katom)
         katom = katom + 1
         go to 10
      else

*        a kekule structure has been generated by including double bonds for 
*        all the atoms. store this structure in kekule(*).

         nkekule = nkekule + 1
         if (nkekule*numkek.gt.memlft) then
            write(iwr,*) ' memory available', memlft, 
     +                   ' memory requied', nkekule*numkek
            call caserr(' memory ran out in subroutine makekekule')
         end if
         do m = 1, numkek
            kekule((nkekule-1)*numkek+m) = klist(m)
         end do
      end if

      return
      end

************************************************************************
      subroutine printkconf(kconf,melec1,nconf,iwr)
      implicit REAL (a-h,o-z)
      dimension kconf(melec1,nconf)
      do i=1,nconf
        write(iwr,*)'in printkconf',i,kconf(melec1,i)
      enddo
      end
      subroutine copi2(ifrom,ito,nelem,incr)
c...
c...  copy integer arrays. the copy-order is reversed by a negative 
c...  increment. (proper choice of ifrom, ito and incr ensures
c...  proper copying if ifrom and ito belong to the same array)
c...  ifrom and ito must be called with their start adresses
c...
      implicit REAL (a-h,o-z), integer (i-n)
c
INCLUDE(../m4/common/iofile)
c
      dimension ifrom(*), ito(*)
c
       write(iwr,'(A,36I3)') 'ifrom::',(ifrom(jj),jj=1,nelem)
       write(iwr,'(A,36I3)') 'ito::',(ito(jj),jj=1,nelem)
      if (incr.lt.0) then
         jfirst = nelem
         jlast  = 1
      else
         jfirst = 1
         jlast  = nelem
      end if
      write(iwr,*)jfirst,jlast,incr
      do 10 icp=jfirst,jlast,incr
         ito(icp) = ifrom(icp)
10    continue
c
       write(iwr,'(A,36I3)') 'ito::',(ito(jj),jj=1,nelem)
      return
      end

