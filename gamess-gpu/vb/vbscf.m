c  $Author: jvl $
c  $Date: 2014-09-15 17:31:14 +0200 (Mon, 15 Sep 2014) $
c  $Locker:  $
c  $Revision: 6296 $
c  $Source: /c/qcg/cvs/psh/GAMESS-UK/vb/vbscf.m,v $
c  $State: Exp $
c  $Log: vbscf.m,v $
c  Revision 1.86  2008-01-31 09:42:54  jvl
c  Super Hybrid fixed (nsa)
c  /marcin
c
c  Revision 1.85  2007/10/18 16:45:34  jvl
c  Problem with clearing nodes on Huygens fixed (inside subroutine vberr).
c  Some cosmetic changes.
c  Hopefully everything works :)
c  /marcin
c
c  Revision 1.82  2007/09/18 10:06:55  jvl
c  VB Parallel bugs fixed for huge VB calculations.
c  HSMAT switched to igmem_alloc. Probably more fixes needed.
c  /marcin
c
c  Revision 1.81  2007/07/04 12:37:56  jvl
c  Introducing new VB switch - VBMO, first step to separate output file for molden.
c  Now, upon request, separate file with VB orbitals (with proper format) can be created,
c  called out.vbmo.<title_from_VBMO_switch>. Check vbin.m comment/GAMESS manual for more.
c  Major bugfixes of dimensions in makepsi0 subroutine:
c   - serial + CORE XX
c   - parallel + bigger basis set than 6-31g + CORE switch
c  Minor changes in igmem_alloc_inf names.
c  /marcin
c
c  Revision 1.80  2007/03/20 14:49:31  jvl
c  Pretty major overhaul of the VB code
c  say 50% is now dynamic using igmem_alloc and, oh wonder,
c  all the examples are still checking out OK, including a new resonating
c  super hybrid. (will follow)
c
c  Revision 1.79  2007/03/15 14:03:09  jvl
c  update super hybrid oprion in vb (ignore ao per atom)
c  added damping
c  clean up
c
c  Revision 1.78  2007/02/20 15:21:16  jvl
c  The array atomao (hybrid per ao) made the hybrid code a lot simpler.
c  Unfortunately now we need the possibility to let an ao belong to two hybrids.
c  Thus packing was in order and a special function (ochatao) was made to handle atomao
c  Further some restrictions were removed; now resomating BLW may be used (:-))
c
c  Revision 1.77  2006/12/13 21:56:08  rwa
c  First check in of the experimental VB response property code.
c
c  Revision 1.76  2006/04/27 19:26:58  jvl
c  Upon finding that the orthogonalisations in VB are inconsistent
c  a cleanup was started, which is no way fnished.
c  natorb is doing all sort of non related stuff and makepsi0 calls natorb
c  which is stupid, also frozen core gradients, which are rather trivial
c  for VB should be properly added, etc. but we have to do some work again
c
c  Revision 1.75  2006/04/18 16:35:34  jvl
c  - added basis option to allow Vb to be run in anothr basis.
c  - cleaned up various bits of VB, like core partitioning
c  - adapted get1e
c
c  Revision 1.74  2006/01/13 17:56:48  jmht
c  Merging in changes from the release branch
c
c  Revision 1.73  2005/09/23 14:34:56  jmht
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
c  Revision 1.72.2.3  2005/10/28 13:43:33  jvl
c  added text to indicate from where a davidson is called
c
c  Revision 1.72.2.2  2005/07/19 06:52:20  wab
c
c  min0/max0 changes (see m4 checkin)
c
c  Revision 1.72.2.1  2005/06/24 05:42:38  wab
c  Second attempt to do this after the bug in guess that caused apparent
c  failures to be diagnosed yesterday has been fixed
c
c  Revision 1.72  2005/04/22 11:07:55  jvl
c  All vb criteria are now in common vbcri and may be set
c  in the input. Also their defaults are now consistent
c
c  Revision 1.71  2005/03/31 15:00:00  jvl
c  - corrected pert option, which was screwed earlier
c  - include file hsinfo
c
c  Revision 1.70  2005/03/25 23:35:47  jvl
c  fixed core partitioning error found in i8 version
c
c  Revision 1.69  2005/03/17 16:28:23  jvl
c  fixed brmax and pert again
c
c  Revision 1.68  2005/03/07 14:08:08  jvl
c  - made sure that hybrids and core work => change small to any
c    (in old small calculations cores would be OK for hybrids (within small)
c    now for bigger calculations this is not the case anymore
c  - tightened critor and orthogonalisation criteria. 10**-10 is not good enough
c    anymore => 10**-14 or 15
c
c  Revision 1.67  2005/02/05 18:04:26  wab
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
c  Revision 1.66  2005/01/19 15:52:49  rwa
c  Minor bug fix memorymanagment for 1-el dens. mat.
c
c  Revision 1.65  2005/01/19 15:46:21  rwa
c  Minor bugfix.  The VB 1-el density matrix was wrong
c  when the wavefunction consists of 1 determinant.
c
c  Revision 1.64  2005/01/05 10:22:53  jvl
c  Bugfix concerning erasing scratch information between consequetive VBSCF
c  runs (for PCM). Further extension of the PCM model with a tweak input
c  option (JJE).
c
c  Revision 1.63  2004/08/24 08:20:34  jvl
c  Dimension bugfix at the end of vbscf (JJE)
c
c  Revision 1.62  2004/08/19 15:52:33  jvl
c  Bugfix in virtual fock space. Bugfix in fock canonicalise option. Better
c  ebi print out. Parameter change (JJE)
c
c  Revision 1.61  2004/05/25 12:27:25  jvl
c  little dimensioning problem still in virtual; n was reset to full size too late (harmonic)
c
c  Revision 1.60  2004/05/24 22:18:21  jvl
c  corrected harmon and cleaned up (ns was used wrongly so canonicalisation failed)
c
c  Revision 1.59  2004/05/19 22:57:45  jvl
c  Brought harmonic to VB; If hamonic is specified the cartesion
c  components are projected out of the occupied and virtual orbitals
c  yielding a harmonic vbscf and vbci (the latter may of may not be what you want)
c
c  Revision 1.58  2004/03/15 07:44:25  jvl
c  bugfix in normalisation of vectors which are specified manually (JJE).
c
c  Revision 1.57  2004/01/29 13:01:39  jvl
c  bugfix in iteration print (JJE)
c
c  Revision 1.56  2004/01/29 12:58:36  jvl
c  Extended TURTLE with real brillouin criterion, which is now default. Fixed
c  some bugs regarding perturbation option. Changed print-out of iterations:
c  brm value added. (JJE)
c
c  Revision 1.55  2004/01/29 11:44:41  jvl
c  changed iteration format
c
c  Revision 1.54  2003/11/09 14:53:41  jvl
c  the layout of ed7 in vb was partly redone every time in an geometry optimisation
c  This cause errors and a possible oveflow of ed7 ....
c
c  Revision 1.53  2003/11/07 16:43:35  jvl
c  fixed two bugs (JJE) one concerning the use of an uninitialised matrix. And
c  the generation of 2-el dens matrix and lagrangian in the case of geometry
c  optimisation.
c
c  Revision 1.52  2003/10/08 17:31:14  jvl
c  several bugfixes (JJE 8/10/2003).
c
c  Revision 1.51  2003/08/28 14:26:51  jvl
c  fixed a bug in array size for igroup (JJE 28/8/2003).
c
c  Revision 1.50  2003/06/10 16:43:53  jvl
c  Addition of option p0core (JJE 10/6/2003).
c
c  Revision 1.49  2003/04/25 13:18:37  jvl
c  Fixed problem with storing vectors on wrong position. Fixed problem with
c  virsig (switched off virsig) and changed several real declarations to
c  REAL declarations (jje).
c
c  Revision 1.48  2003/04/14 05:50:34  jvl
c  Changed output of vbscf, so that it resembles the original output. Natural
c  orbitals arfe no longer printed by default, only the occupation numbers
c  are printed instead (jje 14/4/2003)
c
c  Revision 1.47  2003/03/26 13:46:05  jvl
c  added array idoub as formal parameter; It was declared but not passed which worked
c  depending on the compiler supporting automatic arrays and depending on luck
c
c  Revision 1.46  2003/03/04 14:16:56  jvl
c  (re)fixed bug on sgi. In subroutine makepsi0 something went wrong with
c  reordering of the orbitals, but only on sgi. Put all instructions in new
c  subroutine swapvec, which resolved the problem (jje).
c
c  Revision 1.45  2003/02/28 16:56:37  jvl
c  The number of determinants was calculated wrongly in buildnat. Fixed (jje)
c
c  Revision 1.44  2003/02/28 10:58:33  jvl
c  Fixed bug in subroutine makepsi0. Lowered vbscf.m deck from -O3 to -O2
c  for sgi systems (JJE).
c
c  Revision 1.43  2003/02/19 17:55:57  jvl
c  fixed bug in subroutine overlon: number of determinants was calculation
c  wrong (jje 19/2/2003)
c
c  Revision 1.42  2003/02/18 17:17:17  jvl
c  A lot of changes, in general:
c  - In vbci the call to bsmat for natural orbitals is no longer necessary
c  - The suboptions pert and fock for optimise in scf are implemented
c  - The subroutine vbscf is completely restructured and is more logical
c  - Addition of the makepsi0, which uses a small extra 2-electron,
c    transformation
c    JJE - 18 feb 2003
c
c  Revision 1.41  2002/11/25 15:25:31  jvl
c  There was a restriction to 500000 words for 4-index. Not very relevant
c  and causing multipass mode, which unfortunately crashes in parallel runs.
c  Eliminated restriction and added error-message for multi-pass in parallel
c
c  Revision 1.40  2002/11/18 17:35:28  jvl
c  Modifications to add servec, in VB
c  putqvb is now just interface to putq, rdistr requires maximum dimension
c  schmids has 'norm' keyword to denote, if a normalisation is required
c  and few other routines made more straightforward.
c
c  Revision 1.39  2002/09/12 15:24:04  jvl
c  corrected a confused use of ed7 in vbdiis
c
c  Revision 1.38  2002/09/11 13:23:46  jvl
c  Fixed a bug only showing up in parallel version (JJE). introduced when
c  common block was replaced by an include statement.
c
c  Revision 1.37  2002/09/06 16:56:50  jvl
c  Added external fget for linux in subroutine fockbuild (JJE).
c
c  Revision 1.36  2002/09/05 14:46:16  jvl
c  Changes made for new opti option. subroutine optcorevb added. It calculates
c  Fock matrix. Is not used yet (JJE).
c
c  Revision 1.35  2002/09/04 12:44:09  jvl
c  typo in ed7 addressing ; corrected (thanks to Jeroen)
c
c  Revision 1.34  2002/08/15 13:20:22  jvl
c  *-1.0d0 is a syntax error according to IBM; fixed
c  outstanding rndlmo_vb fix
c
c  Revision 1.33  2002/05/28 15:07:49  jvl
c  General cleanup
c  got rid of all fortran-files (now part of ed7 is used)
c  handle 3 more common's through include
c  make vector sections in vb more logical
c  added print parameter to getqvb
c  ....
c
c  Revision 1.32  2002/05/21 20:44:31  jvl
c  go a lot more relaxed about near singular metric; if > 1.0d-20 only warn
c
c  Revision 1.31  2002/05/03 12:26:02  jvl
c  Cleaned up getqvb/putqvb nouw adapted vectors are ok as inpu
c  added extra vectors combine option
c  made #its in virtual localisation input parameter
c  poporb does not compile well with high optimisation on sgi
c
c  Revision 1.30  2002/04/28 22:30:29  jvl
c  Added virtual options for vb-orbital optimisation
c  in the process an adapted version of poporb is made
c
c  Revision 1.29  2002/04/25 17:21:34  jvl
c  added checks on negative norm on vector in davidson (= overcomplete sywstem)
c  changed a few (print)flags
c
c  Revision 1.28  2002/03/06 22:33:24  jvl
c  some fortran errors flagfged by hp + osinv clash vb -mopac resolved
c
c  Revision 1.27  2002/02/10 21:07:01  jvl
c  cleaned up printing and made comparison ci/scf possible
c
c  Revision 1.26  2002/02/07 12:33:32  jvl
c  fixed bugs + added print of orthog h-matrix
c
c  Revision 1.25  2002/02/01 12:39:31  jvl
c  Bugfix: nactiv was added to iocvi-array, while function mv1v2 expected
c  iocvi without nactiv (number of active orbitals). JJE.
c
c  Revision 1.24  2001/12/17 17:28:50  jvl
c  few extra bugfixes found when porting to linuxppc
c
c  Revision 1.23  2001/11/22 15:56:19  jvl
c  corrected determinant wiping (did not work for branching diagrams)
c  and changed isign usage and arrays which were dynamic
c
c  Revision 1.22  2001/11/13 11:23:25  jvl
c  fixed reduction of double occurring orbitals in redorb
c
c  Revision 1.21  2001/11/05 16:35:55  jvl
c  Two fixes (JJE) :
c  a. The program no longer generates a fatal error during the first iteration
c     on a phase error when DIIS is used. It generates a warning if iteration=1
c     and an error when iteration>1.
c  b. After the creation of the virtual orbitals the subroutine redorb removes
c     the redundant orbitals. In some cases the orbitals were removed twice because
c     the checking procedure wasn't correct (it only decreased the redundant
c     orbital counter by 1).
c
c  Revision 1.20  2001/10/24 21:45:40  jvl
c  corrected few problems in geometry optimisation with VB and added initial stuff
c  to limit output later..
c
c  Revision 1.19  2001/10/19 16:44:07  jvl
c  forgot return/end when removing subr. maporb ; fixed
c
c  Revision 1.18  2001/10/19 16:30:23  jvl
c  clean hybrid (default) now tries to make hybrid exact
c  root broadcasts orbitals + excitation info for parallel jobs
c  hybrids are checked for consistency
c  infato now contains fragment/ao (atomao) numbers making lots of routines easier
c  some routines had name-changes (clao=> clvecvb)
c  eliminated confusing excitation prints
c
c  Revision 1.17  2001/10/15 15:32:42  jvl
c  major parallel overhaul; makes asynchronous switchable
c  no luck yet however ...........
c
c  Revision 1.16  2001/09/19 16:02:45  jvl
c  some corrections resulting from DEC compile
c
c  Revision 1.15  2001/09/13 21:31:17  jvl
c  since it proved impossible to use the message length as eof indicator, now mword=-1 serves this purpose
c  this requires a few changes .....
c
c  Revision 1.14  2001/09/07 14:28:12  jvl
c  common stats in vb clashed with mpi; => stats_vb
c  so about 10 minutes after installing gamess on sara, the first bugfix
c
c  Revision 1.13  2001/07/06 15:20:58  jvl
c  corection of print improvement
c
c  Revision 1.12  2001/07/04 21:21:13  jvl
c  various cosmetic improvements; extended super hybrid option
c  p.s. mchines seems to think all files have changed ...
c
c  Revision 1.11  2001/07/03 14:13:18  jvl
c  fixed a few bugs in parallel code and allowed non-usage of peigss in vb
c  (i.e. peigss specifies if peigas is to be used in vb)
c
c  Revision 1.10  2001/07/02 09:16:05  jvl
c  added ignore/exclude to super hybrid
c  added frozen core contributions properly
c  bit of printing fixes
c
c  Revision 1.9  2001/06/30 13:15:48  jvl
c  made vb input bit clearer + bugfixes
c
c  Revision 1.8  2001/06/27 15:28:48  jvl
c  Changed vector printing
c  fixed various bugs (related to getin2)
c  added frozen en frozen hybrids options (experimental)
c
c  Revision 1.7  2001/06/21 16:38:21  jvl
c  added super hybrid option (experimental and perhaps useless) and a bit of cleaning up
c
c  Revision 1.6  2001/06/12 12:21:19  jvl
c  - Removed several bugs
c  - Improved parallel 4-index transformation including changes in getin2
c  - Added timing for VB routines
c
c  Revision 1.5  2001/06/04 21:59:20  jvl
c  added vectors colmbine option + cleaned vb print
c
c  Revision 1.4  2001/04/14 19:40:09  jvl
c
c  changes b from linuxppc port + some vb fixes
c
c  Revision 1.3  2001/02/08 14:34:37  jvl
c  Adapted parallelisation for GAMESS-UK.
c  An MPI and GA version have been created
c
c  Revision 1.2  2000/03/31 14:18:27  jvl
c  Changed vnorm to vnrm
c  Changed get1e to use getmat
c  Changed scopy to icopy
c  Improved on information about titles
c  Fokke Dijkstra
c
c  Revision 1.1  2000/02/16 12:20:46  jvl
c  adding vb files to gamess repository
c
c Revision 1.10  1998/02/16  14:22:44  fokke
c added counter process to parallel implementation
c added output of energies and interactions of single determinants
c
c Revision 1.9  1997/07/01  13:55:19  fokke
c added save in schmidt (trouble on SGI 7.2 compiler)
c
c Revision 1.7  1997/05/22  12:53:09  joop
c changed imove to icopy
c
c Revision 1.7  1997/05/22  12:53:09  joop
c changed imove to icopy
c
c Revision 1.6  1997/05/22  11:31:48  joop
c Added include macro.cpp
c
c Revision 1.5  1997/01/02  16:21:22  joop
c MPI changes for SGI + 1 bugfix
c
c Revision 1.4  1996/11/07  17:10:01  joop
c SiGr => SGI + bugfixje
c
c Revision 1.3  1996/10/29  15:32:07  joop
c rcs info + parallel + labeled matrix elements
c
c
      subroutine vbscf(v,maxvec,ci)
c
      implicit REAL (a-h,o-z), integer (i-n)
c
c...  control routine for orbital optimisation
c...  v  : orbitals
c...  ci : ci vector passed from parent subroutine
c
INCLUDE(../m4/common/sizes)
INCLUDE(../m4/common/iofile)
INCLUDE(common/turtleparam)
      dimension v(*),ci(*)
c
INCLUDE(common/vbdiist)
INCLUDE(common/tractlt)
      common /bypass/ index4,ihmat,idavid
INCLUDE(common/hsinfo)
INCLUDE(common/vbcri)
c
c...  for dynamic memory allocation debugging purpose, contains IGMEM_ vars.
INCLUDE(../m4/common/gmemdata)
INCLUDE(../m4/common/gmempara) 
INCLUDE(../m4/common/vcore) 
c
INCLUDE(common/brill)
c
INCLUDE(common/scftvb)
INCLUDE(common/ffile)
INCLUDE(common/splice)
c
      common /davcom/ eshift,thresj,thresd,maxcyc,maxdav,maxsel,iprind,
     &                alter,firsts
      logical alter,firsts
c
INCLUDE(common/infato)
INCLUDE(common/vbequiv)
c
INCLUDE(common/vbpert)
c
INCLUDE(../m4/common/mapper)
INCLUDE(../m4/common/restri)
INCLUDE(../m4/common/timeperiods)
c
INCLUDE(common/twice)
INCLUDE(common/first_vb)
INCLUDE(common/hsmattype)
c
INCLUDE(common/aivb)
INCLUDE(common/vbtimess)
c
      common /txema/ struc,ivec
      logical struc
c
      common/checkblock/ icheckblock
      common/gjs   /nirr,mult(8,8),isymao(maxorb),isymmo(maxorb)
     1              ,irr,iss,nstart(8),nfin(8),nmc
c
      character*10 zstr
      logical oroot
      integer idumpr(5,2)
c
c...  ncol2 is # of vectors before redundant ones are deleted
c
      save ncol2
c
      icheckblock = 0
      call init_vbdiis
      swave = 0.0d0                    
      evbpre = 0.0d0
c
      ntscf = ntscf + 1
c
      tvirtu = 0.0d0
      t4indx0 = 0.0d0
      n4indx0 = 0
      t4indx = 0.0d0
      n4indx = 0
      tmatre0 = 0.0d0
      nmatre0 = 0
      tmatre = 0.0d0
      nmatre = 0
      tdavid = 0.0d0
      tlagr  = 0.0d0
      tgtr   = 0.0d0
      minfock = 999999999
      maxfock = 0
      minfockd = 999999999
      maxfockd = 0
      minpert = 999999999
      maxpert = 0
      minvar = 999999999
      maxvarp = 0
c
      kmaxmem = igmem_max_memory()
      kmemscr = kmaxmem/10
      kmemleft = kmaxmem - kmemscr
      ks   = igmem_alloc_inf(nbasis*(nbasis+1)/2,'vbscf.m','vbscf', 
     &                       'ks',IGMEM_DEBUG)
      kscr = igmem_alloc_inf(kmemscr,'vbscf.m','vbscf','kscr',
     &                           IGMEM_DEBUG)
      call get1e(Q(ks),dummy,'s',Q(kscr))
      if (ofirst.and.iprins.gt.50) then
         if (iprins.gt.10000) then
            write(iwr,1000)
1000        format(/,' the ao metric :')
            call tripri(Q(ks),nbasis)
         end if
         call sigorb(Q(ks),nbasis)
      end if
c
c...  nitscf = iteration # in the scf process
c...  scfconv= true if the optimisation has converged
c...  nactiv = # active, i.e. occupied, mo's
c...  ncol   = # orbitals (core,active,virtual), without redundant
c...  ncol2  = # orbitals (core,active,virtual), with redundant
c...             (redorb will reduce ncol)
c...  nsa    = # orbitals in 4-index (active+virtual)
c
      call getqvb(v,nbasis,ncol,isecv,'nopr')
      call normt(v,nbasis,ncol,Q(ks),Q(kscr))
c...     make sure equivalence restrictions are obeyed 
      if (equiv) call primequiv(v,nbasis,1,'occ','at start')
c
      call putqvb(v,nbasis,ncol)
      call vbfirst(v,maxvec)
      call gmem_free_inf(kscr,'vbscf.m','vbscf','kscr')
      call gmem_free_inf(ks,'vbscf.m','vbscf','ks')
c
      nitscf = 0
      scfconv = .false.
      npert = 0
      evb = 0.0d0
c
      nactiv = nsa
c...  (this nsa was the # active orbitals in the active statement)
      nsa = ncol - ncore
      nvirt = nsa - nactiv
      ncol2 = ncol
c
c...  start of iterative loop
c
1     nitscf = nitscf + 1
c
c...   perform allowed orthogonalisations (if ortscf)
c
      mm = ncore + nactiv + nvirt
c      call prvc14(v,mm,nbasis,v,'o','l')
      if (ortscf) 
     &   call vbcanon(v,nbasis,ncore+nactiv,100,'occ',Q(1))
*
c...     make sure equivalence restrictions are obeyed 
c
      if (equiv) call primequiv(v,nbasis,1,'occ','in scf iteration')
c
c...  calculate psi(0)
c       
_IF(debug)
        write(iwr,*) '(1) vbscf'
_ENDIF
*
        evbpre = evb
*
10      call makepsi0(v,ci)
*     
        if (ovbdiis) then
           if (evbpre.lt.evb) then
*     
              write(iwr,20)     
20            format(100('-'))
              write(iwr,30) evbpre,evb
30            format(6x,'==> DIIS failed  previous E = ',f22.14,
     +              ' last E = ',f22.14,' <==')
              write(iwr,20)
              kvb7_vnew = kscra7vb('kvb7_vnew',nbasis*nactiv,'r','r')
              call rdedx(v(ncore*nbasis+1),nbasis*nactiv,kvb7_vnew,num8)
              ovbdiis = .false.
              call init_vbdiis
              go to 10
*      
           end if
        end if
c
        if (fckp) then
           a1 = cpulft(1)
           call fockvb(Q(1))
           tfock = tfock + cpulft(1) - a1
           nfock = nfock + 1
        endif
_IF(debug)
        write(iwr,*) '(2) vbscf'
_ENDIF
c
c...  perform orthogonalisations as defined in forbex
c
        call restror(v(ncore*nbasis+1),nsa,nbasis,Q(1))
*
_IF(debug)
        write(iwr,*) '(3) vbscf'
_ENDIF
c
c...  construct virtual space (see virtual)
c
        a1 = cpulft(1)
        call start_time_period(TP_VB_VIRT)
c
c...  use special techniques to get a properly shaped virtual
c...  space (see "virtual")
c
*
*      why can't I see the reason of allocating ksao here (zahid)
*
c        ksao = igmem_alloc_inf(nbasis*(nbasis+1)/2,'vbscf.m','vbscf',
c     &                         'ksao',IGMEM_DEBUG)
_IF(debug)
        write(iwr,*) '(4) vbscf'
_ENDIF
*
        kocc = 1
        maxvirt = maxvec - (nscf+ncore)
        call virtual(v(kocc),v(kocc + (nscf+ncore)*nbasis),
     &               nbasis,ncore,ncore+nscf,nvirt,maxvirt)
        if (super_hybrid) nsa = ncol - ncore
c        call gmem_free_inf(ksao,'vbscf.m','vbscf','ksao')
c
c...  remove redundant virtuals
c
c...  the vectors can be changed by redorb ! (redundant
c...  virtuals are thrown away, iex is adapted) nredund
c...  sits in /brill/. ncol is saved in ncol2
c
      if (reduce) then
        if (equiv.and.nitscf.le.0) print *,' Remco : redorb is fatal'
        if (nactiv.ne.nsingly+ndoubly) 
     1     call caserr('nactiv ne nsingly+ndoubly')
_IF(debug)
        write(iwr,*) '(5) vbscf - reduce'
_ENDIF
        ksao = igmem_alloc_inf(nbasis*(nbasis+1)/2,'vbscf.m','vbscf',
     &                         'ksao',IGMEM_DEBUG)
        kvcopy = igmem_alloc_inf(nbasis*(nbasis+nactiv+ncore),
     &                           'vbscf.m','vbscf','kvcopy',
     &                           IGMEM_DEBUG)
        kiocvi = igmem_alloc_inf((nactiv+ncore-1)/nipw()+1,'vbscf.m',
     &                           'vbscf','kiocvi',IGMEM_DEBUG)
        kiset = igmem_alloc_inf(nbasis,'vbscf.m','vbscf','kiset',
     &                          IGMEM_DEBUG)
        kiact = igmem_alloc_inf(nactiv,'vbscf.m','vbscf','kiact',
     &                          IGMEM_DEBUG)
        kidoc = igmem_alloc_inf(nactiv,'vbscf.m','vbscf','kidoc',
     &                          IGMEM_DEBUG)
        kisoc = igmem_alloc_inf(nactiv,'vbscf.m','vbscf','kisoc',
     &                          IGMEM_DEBUG)
        kivir = igmem_alloc_inf(nbasis,'vbscf.m','vbscf','kivir',
     &                           IGMEM_DEBUG)
        nam = nbasis + nscf
        khmoao = igmem_alloc_inf(nam*(nam+1)/2,'vbscf.m','vbscf',
     &                           'khmoao',IGMEM_DEBUG)
        kiex2 = igmem_alloc_inf(nactiv*(nbasis+1),'vbscf.m','vbscf',
     &                          'kiex2',IGMEM_DEBUG)
_IF(debug)
        write(iwr,*) '(55) vbscf - reduce'
_ENDIF
*
*       strangely enough the code was working even without getting
*       the s matrix (zahid) 
*
        call get1e(Q(ksao),dumm,'s',Q(ksao))
*
        call redorb(Q(ksao),v(ncore*nbasis+1),Q(kvcopy),v,
     +              Q(kiocvi),Q(kiset),Q(kiact),
     +              Q(kidoc),Q(kisoc),Q(kivir),
     +              Q(khmoao),Q(kiex2),nbasis,Q(1))
        ncol = ncol2 - nredund
        nsa = ncol - ncore
        nvirt = nsa - nactiv
        call gmem_free_inf(kiex2,'vbscf.m','vbscf','kiex2')
        call gmem_free_inf(khmoao,'vbscf.m','vbscf','khmoao')
        call gmem_free_inf(kivir,'vbscf.m','vbscf','kivir')
        call gmem_free_inf(kisoc,'vbscf.m','vbscf','kisoc')
        call gmem_free_inf(kidoc,'vbscf.m','vbscf','kidoc')
        call gmem_free_inf(kiact,'vbscf.m','vbscf','kiact')
        call gmem_free_inf(kiset,'vbscf.m','vbscf','kiset')
        call gmem_free_inf(kiocvi,'vbscf.m','vbscf','kiocvi')
        call gmem_free_inf(kvcopy,'vbscf.m','vbscf','kvcopy')
        call gmem_free_inf(ksao,'vbscf.m','vbscf','ksao')        
_IF(debug)
        write(iwr,*) '(555) vbscf - reduce'
_ENDIF
      end if
_IF(debug)
        write(iwr,*) '(6) vbscf'
_ENDIF
c
c...  dump vectors (including virtuals) for orthopt
c
        nn = ncore + nactiv + nvirt
        kvb7_vv = kscra7vb('kvb7_vv',nn*nbasis,'r','w')
        call wrt3(v,nn*nbasis,kvb7_vv,num8)
c
c...  print the orbitals / excitation patterns the first iter.
c
      if (nitscf.eq.1) then
         nn = nactiv
         if (iprins.gt.1) nn = ncore + nactiv
         if (iprins.gt.10) nn = ncore + nactiv + nvirt
      else
         nn = 0
         if (iprins.ge.5) nn = nactiv
         if (iprins.ge.50) nn = ncore + nactiv + nvirt
      end if
      if (nn.eq.nactiv) then
         kvv = ncore*nbasis + 1
         zstr = 'excl. core'
         kk = 0
      else
         kvv = 1
         zstr = 'incl. core'
         kk = ncore
      end if
c
_IF(debug)
        write(iwr,*) '(7) vbscf'
_ENDIF
      if (nn.ge.nactiv) then
         write(iwr,609) ncore,'frozen','fzc',(i,i=1,ncore)
         write(iwr,609) ndoubly,'doubly','doc',
     1                 (idoubly(i)+ncore,i=1,ndoubly)
         write(iwr,609) nsingly,'singly','voc',
     1                 (isingly(i)+ncore,i=1,nsingly)
         write(iwr,610) nvirt,nbasis
609      format(1x,'+++++',i6,1x,a6,' occupieds (',a3,'):',(t36,15i4))
610      format(1x,'+++++',i6,1x,' virtuals / ',i6,' basis functions')
         if (nsingly+ndoubly.ne.nactiv) 
     l       call caserr('nsingly-ndoubly - nactiv  clash')
c
         write(iwr,601) zstr
601      format(1x,/,'         ============================= ',
     1             /,'          real excitations ',a10,
     2             /,'         ============================= ')
         do i=1,nscf
           if (nex(i).gt.0) write(iwr,602) (iex(j,i)+kk,j=1,nex(i)+1)
602        format('       mo',i3,' ==>',15i4,/,(16x,15i4))
         end do
c
c... check
c
         kk = 0
         do i=1,nscf
            do j=1,nex(i)+1
               do k=j+1,nex(i)+1
                  if (iex(j,i).eq.iex(k,i)) then
                     kk = kk + 1
                  end if
               end do
            end do
         end do
         if (kk.gt.0) then
            write(iwr,'(a,i5)') ' *** found double arbitals # ',kk
            call caserr(' program error - double excitations')
         end if
c
         if (super_hybrid.and.nn.le.ncore+nactiv.and.igno_sel.ne.99)then
            write(iwr,606)
606         format(1x,/,'         ============================= ',
     1                /,'          super_hybrid - virt vs. ao ',
     2                /,'         ============================= ')
            kkkv = (ncore+nactiv)*nbasis + 1
            kkk = 1
            do i=1,nvirt
               call abmax(v(kkkv),nbasis,1,'sq',1,am,im)
               if (am.ne.1.0d0) call caserr('super_hybrid funny')
               idumpr(kkk,1) = i+nactiv+kk
               idumpr(kkk,2) = im
               if (kkk.eq.5) then
                  write(iwr,607) (idumpr(k,1),idumpr(k,2),k=1,kkk)
                  kkk = 0
               end if
               kkk = kkk + 1
               kkkv = kkkv + nbasis
            end do
            kkk = kkk - 1
            if(kkk.ne.0)write(iwr,607) (idumpr(k,1),idumpr(k,2),k=1,kkk)
607         format(7x,10(:' (',i4,' - ',i4,')  '))
         end if
c
         if (nn.eq.ncore+nactiv) write(iwr,603) ncore
         if (nn.lt.ncore+nactiv) write(iwr,604)
         if (nn.gt.ncore+nactiv) write(iwr,605) ncore,nvirt
603      format(1x,/,'             ===================== ',
     1             /,'                  real vectors',
     1             /'                 incl',i3,' core',
     2             /,'             ===================== ')
604      format(1x,/,'             ===================== ',
     1             /,'                  real vectors',
     2             /,'             ===================== ')
605      format(1x,/,'             ===================== ',
     1             /,'                  real vectors',
     1             /'            incl',i3,' core and',i3,' virt',
     2             /,'             ===================== ')
         call prvc(v(kvv),nn,nbasis,v,'o','l')
      end if
c
c...  check hybrids
c

_IF(debug)
        write(iwr,*) '(8) vbscf'
_ENDIF

      if (hybry) call clvecvb(v(ncore*nbasis+1),nactiv+nvirt,nbasis,
     1                        'check','null','in vbscf')
_IF(parallel)
c
c...  to be sure broadcast orbitals en info; to get all noses aligned
c
      call pg_brdcst(8,v,8*nbasis*(ncore+nactiv+nvirt),0)
      nn = maxact*(1+maxex+1+maxex+maxex+1)+1+mxorbvb+1+
     1     (maxact+1)*maxact+maxact+7+maxact
      nn = nn*8/nipw()
      call pg_brdcst(9,nex,nn,0)
_ENDIF
c
c...  save orbitals to disk as these are the vbscf orbitals
c
        call putqvb(v,nbasis,nactiv+ncore)
c
        call flushbuffer
_IF(peigss)
        call cleargas()
_ENDIF
c...  the stuff before had to be done for consistent orthog etc.
        if (scfconv) goto 9999
c
        a1 = cpulft(1) - a1
        tvirtu = tvirtu + a1
        call end_time_period(TP_VB_VIRT)
c
c...  transform the integrals (2nd transformation including virtuals)
c
        a1 = cpulft(1)
        call start_time_period(TP_VB_TRAN)
        nmc = nsa + ncore 
        if (scfconv) nsa = nactiv
        lenact = nsa*(nsa+1)/2
c
        ks = igmem_alloc_inf(lenact,'vbscf.m','vbscf','ks',IGMEM_DEBUG)
        kh = igmem_alloc_inf(lenact,'vbscf.m','vbscf','kh',IGMEM_DEBUG)
        ng = lenact*(lenact+1)/2+1
_IF(debug)
        write(iwr,*) '(9) vbscf'
_ENDIF
c
c...  (first integral = (00/00) = 0.0 => + 1)
c...  don't use to much core in tran as virtual io may interfere
c...  with normal io; If the user wants to commit suicide, let him (jvl,2002)
c
        call transformvb(Q(ks),Q(kh),v)
_IF(debug)
        write(iwr,*) '(10) vbscf'
_ENDIF
        if ( nitscf.eq.1 .and. n2int+1 .ne. lenact*(lenact+1)/2+1) then
          write(iwr,'(a,I10,a,I10,a)')
     &    ' Truncating nr of 2-el integrals (Brillouin) to ',
     &    n2int+1,' from ',lenact*(lenact+1)/2 + 1,'.'
        end if

        kg = igmem_alloc_inf(n2int+1,'vbscf.m','vbscf',
     &                       'kg',IGMEM_DEBUG)
        call getin2(Q(kg))
        call clredx
c
c... print integrals over orbitals if requested
c   
       if (iprint.ge.50) then
          write(iwr,*) ' 1-electron integrals over orbitals'
          call tripri(Q(kh),nsa)
          write(iwr,*) ' overlap matrix between orbitals'
          call tripri(Q(ks),nsa)
          if (iprint.ge.100000) then
            write(iwr,*) ' 2-electron integrals'
            call tripri(Q(kg+1),lenact)
          end if
       end if
c
        a1 = cpulft(1) - a1
        t4indx = t4indx + a1
        n4indx = n4indx + 1
        call end_time_period(TP_VB_TRAN)
        call flushbuffer()
c
c...  calculate brillouin-structure interaction matrix
c
        a1 = cpulft(1)
        call start_time_period(TP_VB_ME)
        call bsmat(Q(ks),Q(kh),Q(kg),
     +             ci,nbasis,nactiv,nvirt)
_IF(debug)
        write(iwr,*) '(11) vbscf'
_ENDIF
        a1 = cpulft(1) - a1
        tmatre = tmatre + a1
        nmatre = nmatre + 1
        call end_time_period(TP_VB_ME)
c
c...  diagonalise bi-matrix (davids takes care of parameter switch)
c
        kmaxmem = igmem_max_memory()
        kmemscr = kmaxmem/10
        kmemleft = kmaxmem - kmemscr
        kbi = igmem_alloc_inf(kmemscr,'vbscf.m','vbscf','kbi',
     &                        IGMEM_DEBUG)
        a1 = cpulft(1)
        call start_time_period(TP_VB_DIAG)

_IF(debug)
        write(iwr,*) '(12) vbscf'
_ENDIF
        call davids(Q(kbi),kmemscr,core,ebi,'bi')
        a1 = cpulft(1) - a1
        tdavid = tdavid + a1
        call end_time_period(TP_VB_DIAG)
c
        ebisci=ebi
        ebi=0.0d0
        if (ovbper) then
_IF(debug)
        write(iwr,*) '(13) vbscf - ovbper'
_ENDIF
c
c...  davids only dealt with internal excitations, guess the
c...  other contributions using perturbation theory. nstruc and
c...  nvarp (# variational treated states) temporarily were
c...  interchanged, restore this.
c
           a1 = cpulft(1) 
           call start_time_period(TP_VB_ME)
           nnn = nstruc
           nstruc = nvarp
           nvarp = nnn
           kciint = igmem_alloc_inf(nvarp,'vbscf.m','vbscf','kciint',
     &                              IGMEM_DEBUG)
           kgrh = igmem_alloc_inf(nstruc,'vbscf.m','vbscf','kgrh',
     &                            IGMEM_DEBUG)
           kgrs = igmem_alloc_inf(nstruc,'vbscf.m','vbscf','kgrs',
     &                            IGMEM_DEBUG)
           kdiagh = igmem_alloc_inf(nstruc,'vbscf.m','vbscf','kdiagh',
     &                              IGMEM_DEBUG)
           kdiags = igmem_alloc_inf(nstruc,'vbscf.m','vbscf','kdiags',
     &                              IGMEM_DEBUG)
           kigr_brill = igmem_alloc_inf((nbrill*5-1)/nipw()+1,
     &                                   'vbscf.m','vbscf',
     &                                   'igr_brill',IGMEM_DEBUG)
c...  move variational vector
           call fmove(Q(kbi),Q(kciint),nvarp)
           call vbper(Q(kciint),Q(kbi),Q(kgrh),
     +                Q(kgrs),Q(kdiagh),Q(kdiags),
     +                Q(kigr_brill))

_IF(debug)
        write(iwr,*) '(133) vbscf - ovbper'
_ENDIF
           a1 = cpulft(1) - a1
           tmatre = tmatre + a1
           nmatre = nmatre + 1
           call end_time_period(TP_VB_ME)

           call gmem_free_inf(kigr_brill,'vbscf.m','vbscf','igr_brill')
           call gmem_free_inf(kdiags,'vbscf.m','vbscf','kdiags')
           call gmem_free_inf(kdiagh,'vbscf.m','vbscf','kdiagh')
           call gmem_free_inf(kgrs,'vbscf.m','vbscf','kgrs')
           call gmem_free_inf(kgrh,'vbscf.m','vbscf','kgrh')
           call gmem_free_inf(kciint,'vbscf.m','vbscf','kciint')
        end if
_IF(debug)
        write(iwr,*) '(14) vbscf'
_ENDIF
        ebi=ebisci+ebi
c
c...  update orbitals and check convergence
c
        knewv = igmem_alloc_inf(nbasis*nactiv,'vbscf.m','vbscf','knewv',
     &                          IGMEM_DEBUG)
        kmaxmem  = igmem_max_memory()
        kmemscr  = kmaxmem/10
        kmemleft = kmaxmem - kmemscr
        kscr     = igmem_alloc_inf(kmemscr,'vbscf.m','vbscf','kscr',
     &                           IGMEM_DEBUG)
c
c...  change orbitals
c
        call orbopt(v,nbasis,nactiv,Q(knewv),Q(kbi),
     +              Q(kscr),kmemscr)
        call gmem_free_inf(kscr,'vbscf.m','vbscf','kscr')
        call gmem_free_inf(knewv,'vbscf.m','vbscf','knewv')
        call gmem_free_inf(kbi,'vbscf.m','vbscf','kbi')
        call gmem_free_inf(kg,'vbscf.m','vbscf','kg')
        call gmem_free_inf(kh,'vbscf.m','vbscf','kh')
        call gmem_free_inf(ks,'vbscf.m','vbscf','ks')
*
        if (unitary) then
           ksao = igmem_alloc_inf(nbasis*(nbasis+1)/2,'vbscf.m','vbscf',
     +                            'ksao',IGMEM_DEBUG)
           kscr = igmem_alloc_inf(nbasis*(nbasis+1)/2,'vbscf.m','vbscf',
     +                          'kscr',IGMEM_DEBUG)
*           
           call get1e(Q(ksao),dummy,'s',Q(kscr))
              
           call gmem_free_inf(kscr,'vbscf.m','vbscf','kscr')
*              
           kscr    = igmem_alloc_inf(nbasis*2,'vbscf.m','vbscf','kscr',
     +                                IGMEM_DEBUG)
*              
           call normvc(v(ncore*nbasis+1),Q(ksao),Q(kscr),nbasis,
     +                 nactiv,cridep)
*              
           call gmem_free_inf(kscr,'vbscf.m','vbscf','kscr')
           call gmem_free_inf(ksao,'vbscf.m','vbscf','ksao')
        end if
*
_IF(debug)
        write(iwr,*) '(15) vbscf'
_ENDIF
c
c...  end of iterative loop. Start on top 
c
      goto 1
c
9999  continue
      icheckblock = 1
c
      ks = igmem_alloc_inf(lenact,'vbscf.m','vbscf','ks',IGMEM_DEBUG)
      kh = igmem_alloc_inf(lenact,'vbscf.m','vbscf','kh',IGMEM_DEBUG)

      nsa = nactiv
c
c...  read info from tape
c
      call wr15vb(nelec,nalfa,nconf,ndets,nstruc,maxdet,maxstr,
     &            nwpack,ncoeff,imax,ndum11,ndoub,'read')
      norb = nsa + nvirt
c
      ksao = igmem_alloc_inf(nbasis*(nbasis+1)/2,'vbscf.m','vbscf',
     &                       'ksao',IGMEM_DEBUG)
      ksao2 = igmem_alloc_inf(nbasis*(nbasis+1)/2,'vbscf.m','vbscf',
     &                        'ksao2',IGMEM_DEBUG)
      kmaxmem  = igmem_max_memory()
      kmemscr  = kmaxmem/10
      kmemleft = kmaxmem - kmemscr
      kscr = igmem_alloc_inf(kmemscr,'vbscf.m','vbscf','kscr',
     &                       IGMEM_DEBUG)
c
c...  calculate s-matrix over orbitals (all but virtual orbitals)
c
      call get1e(Q(ksao),dummy,'s',Q(ksao2))
      call fmos(Q(ks),Q(ksao),v(ncore*nbasis+1),
     +          Q(kscr),nsa,nbasis,crilow)

_IF(debug)
        write(iwr,*) '(16) vbscf'
_ENDIF
      call gmem_free_inf(kscr,'vbscf.m','vbscf','kscr')
      call gmem_free_inf(ksao2,'vbscf.m','vbscf','ksao2')
      call gmem_free_inf(ksao,'vbscf.m','vbscf','ksao')
c
c...  keep certain parameters
c
      hsmatt = 'vbfunction'
      oldncol = ncol
      ncol = ncore + nscf
      oldncore = ncore
      ncore = ncore + ndoub
      oldnsa = nsa
      nsa = ncol - ncore
      lenact = nsa*(nsa+1)/2
      oldncurt = ncurt
      ncurt = -1
      struc=.true.
      kmaxmem  = igmem_max_memory()
      kmemscr  = kmaxmem/10
      kmemleft = kmaxmem - kmemscr
      kscr = igmem_alloc_inf(kmemscr,'vbscf.m','vbscf','kscr',
     &                       IGMEM_DEBUG)
*
      call scfina(v,ci,Q(ks),Q(kscr),Q(1))
*
_IF(debug)
        write(iwr,*) '(17) vbscf'
_ENDIF

      call gmem_free_inf(kscr,'vbscf.m','vbscf','kscr')
c
c...  find the orthogonality classes in the integrals
c
c...  why do we even do this??
c
      kmaxmem  = igmem_max_memory()
      kmemscr  = kmaxmem/10
      kmemleft = kmaxmem - kmemscr
      kiort = igmem_alloc_inf(norb,'vbscf.m','vbscf','kiort',
     &                       IGMEM_DEBUG)
      call sym1(Q(ks),norb,Q(kiort),northo,'noprint')

_IF(debug)
        write(iwr,*) '(18) vbscf'
_ENDIF

      kigro = igmem_alloc_inf(5*nconf,'vbscf.m','vbscf','kigro',
     &                        IGMEM_DEBUG)
      kpacd = igmem_alloc_inf(nwpack,'vbscf.m','vbscf','kpacd',
     &                        IGMEM_DEBUG)
      kndet = igmem_alloc_inf(nstruc,'vbscf.m','vbscf','kndet',
     &                        IGMEM_DEBUG)
      kidps = igmem_alloc_inf(ndets*nstruc,'vbscf.m','vbscf','kidps',
     &                        IGMEM_DEBUG)
      kcoef = igmem_alloc_inf(ncoeff*ndets,'vbscf.m','vbscf','kcoef',
     &                        IGMEM_DEBUG)
      kidet = igmem_alloc_inf(nelec*ndets,'vbscf.m','vbscf','kidet',
     &                        IGMEM_DEBUG)
      kjdet = igmem_alloc_inf(nelec*ndets,'vbscf.m','vbscf','kjdet',
     &                        IGMEM_DEBUG)
      kmaxmem  = igmem_max_memory()
      kmemscr  = kmaxmem/10
      kmemleft = kmaxmem - kmemscr
      kscr  = igmem_alloc_inf(kmemscr,'vbscf.m','vbscf','kscr',
     &                        IGMEM_DEBUG)
c
c...  read information from datatape
c
      call inform(nelec,nalfa,ndets,nstruc,Q(kpacd),
     &            Q(kndet),Q(kidps),Q(kcoef),
     &            ncoeff,Q(kigro),ngroup,Q(kiort),
     &            northo,Q(kscr),kmemscr)
_IF(debug)
        write(iwr,*) '(19) vbscf'
_ENDIF
c
c...  generate psi0 on determinant basis
c
      call psi0det(ci,nstruc,Q(kigro),ngroup,
     &             Q(kcoef),ncoeff,ndets,Q(kndet),
     &             Q(kidps),Q(kidet),Q(kjdet),
     &             Q(kpacd),nelec,nalfa,ndettot,
     &             Q(kscr),kmemscr,
     &             'dontsaveondisk')
_IF(debug)
        write(iwr,*) '(20) vbscf'
_ENDIF
c
c...  move original parameters back
c
      struc=.false.
      hsmatt = 'full      '
      ncore = oldncore
      ncol = oldncol
      nsa = oldnsa
      lenact = nsa*(nsa+1)/2
      ncurt = oldncurt
_IF(debug)
        write(iwr,*) '(21) vbscf'
_ENDIF
      call gmem_free_inf(kscr,'vbscf.m','vbscf','kscr')
      call gmem_free_inf(kjdet,'vbscf.m','vbscf','kjdet')
      call gmem_free_inf(kidet,'vbscf.m','vbscf','kidet')
      call gmem_free_inf(kcoef,'vbscf.m','vbscf','kcoef')
      call gmem_free_inf(kidps,'vbscf.m','vbscf','kidps')
      call gmem_free_inf(kndet,'vbscf.m','vbscf','kndet')
      call gmem_free_inf(kpacd,'vbscf.m','vbscf','kpacd')
      call gmem_free_inf(kigro,'vbscf.m','vbscf','kigro')
      call gmem_free_inf(kiort,'vbscf.m','vbscf','kiort')
      call gmem_free_inf(kh,'vbscf.m','vbscf','kh')
      call gmem_free_inf(ks,'vbscf.m','vbscf','ks')
_IF(debug)
        write(iwr,*) '(KONIEC) vbscf'
_ENDIF
      return
      end
c***********************************************************************
      subroutine bcreat(iref,nelec,nalfa,ibrill,nbrill,ifrom,ito)
c
      implicit REAL  (a-h,o-z) , integer   (i-n)
c
      dimension iref(*),ibrill(*)
      ib = 0
      nbrill = 0
c
c...  alpha loop
c
      do 40 j=1,nalfa
         if (iref(j).eq.ifrom) then
            do 5 k=1,nalfa
5           if (ito.eq.iref(k)) goto 40
            nbrill = nbrill + 1
            do 10 k=1,nelec
               ibrill(ib+k) = iref(k)
10          continue
            ibrill(ib+j) = ito
            ib = ib + nelec
         end if
40    continue
c
c...  beta loop
c
      do 80 j=nalfa+1,nelec
         if (iref(j).eq.ifrom) then
            do 7 k=nalfa+1,nelec
7           if (ito.eq.iref(k)) goto 80
            nbrill = nbrill + 1
            do 50 k=1,nelec
               ibrill(ib+k) = iref(k)
50          continue
            ibrill(ib+j) = ito
            ib = ib + nelec
         end if
80    continue
      return
      end

c***********************************************************************
      subroutine bmain(pacdet,idet  ,jdet  ,detcomb,dettot,
     &                 icp ,jcp ,idetps,coeff  ,ig    ,
     &                 igroup,hamil  ,overl ,
     &                 supers,superh,superg,ipos   ,weight,
     &                 g     ,iortho,s     ,scr1   ,scr2  ,
     &                 scr3  ,ndetps,
     &                 ciref,q,lword,nvirt)
 
      implicit REAL  (a-h,o-z) , integer   (i-n)
c
c...  main routine for the brillouin scf procedure
c...  array-lengths are defined in "bsmat".
c
INCLUDE(../m4/common/sizes)
INCLUDE(common/turtleparam) 
INCLUDE(common/c8_16vb) 
INCLUDE(common/hsinfo)
INCLUDE(common/brill)
INCLUDE(common/scftvb)
INCLUDE(common/vbpert)
INCLUDE(../m4/common/iofile)
INCLUDE(common/splice)
INCLUDE(common/tractlt)
_IF(peigss)
INCLUDE(common/gainfo)
      logical oga,pg_create
      parameter (ichunk=10,jchunk=10)
_ENDIF     
      dimension pacdet(*),idet  (*),jdet  (*),detcomb(*),dettot(*),
     &          icp (*),jcp (*),idetps(*),coeff  (*),ig    (*),
     &        igroup(5,*),hamil (*),overl (*),
     &          supers(*),superh(*),superg(*),ipos   (*),weight(*),
     &          g     (*),iortho(*),s     (*),scr1   (*),scr2  (*),
     &          scr3  (*),ndetps(*),q(lword),
     &          ciref(*)
      norb = nsa 
c
c...  find the orthogonality classes in the integrals
c
      call sym1(supers,norb,iortho,northo,'noprint')
c
c...  read information from datatape
c
      call inform(nelec,nalfa,ndets,nstruc,pacdet,ndetps,idetps,coeff,
     &             ncoeff,igroup,ngroup,iortho,northo,q,lword)
c
c...  generate psi0 on determinant basis 
c
      call psi0det(ciref,nstruc,igroup,ngroup,coeff,ncoeff,
     &   ndets,ndetps,idetps,idet,jdet,pacdet,nelec,nalfa,ndtot,q,lword,
     &   'dontsaveondisk')
c
c... check overlap of old wavefunction with new one
c
      call overlon(nelec,nalfa,igroup,ngroup,idet,q,lword)
c
c...  formulate the brillouin wave-function
c
      call brigen(nelec,nalfa,igroup,ngroup,pacdet,coeff,
     &            coeff(ncoeff+1),ciref,ndetps,idetps,idet,jdet,
     &            iortho,northo,q,lword,superh,supers,superg)
c
c
_IF(peigss)
      if (ovbper) then
        n = nvarp
      else
        n = ngroup
      endif
      call cleargas()
      oga = pg_create(8130,n,n,'hmat_ga',ichunk,jchunk,iga_h)
      if (.not.oga) call ga_error
     &                  ('GA: failed to create hmat_ga',1)
      call pg_zero(iga_h)
      oga = pg_create(19130,n,n,'smat_ga',ichunk,jchunk,iga_s)
      if (.not.oga) call ga_error
     &                  ('GA: failed to create smat_ga',1)
      call pg_zero(iga_s)
_ENDIF
c
c
c...  If fock matrix elements will be used in Brillouin state matrix, transform
c...  the fock matrix on ao basis to mo basis and store it on ed7
c
      if (fckp) then
        timing = cpulft(1)
        nn = ncore + nsa
        lennn=nn*(nn+1)/2
        lennb=nbasis*(nbasis+1)/2
        kvec = 1
        kfmo = kvec + nn * nbasis
        kfao = kfmo + lennn
        kdag = kfao + lennb
        kscr = kdag + nn * nbasis
        kvb7_vv = kscra7vb('kvb7_vv',nn*nbasis,'r','n')
        kvb7_fcao = kscra7vb('kvb7_fcao',lennb,'r','n')
        if ((kvb7_vv.lt.0).or.(kvb7_fcao.lt.0)) then
          call caserr('no vectors/fock matrix on ed7 in bmain')
        else
          kvb7_vv = kscra7vb('kvb7_vv',nn*nbasis,'r','r')
          call rdedx(q(kvec),nbasis*nn,kvb7_vv,num8)
          kvb7_fcao = kscra7vb('kvb7_fcao',lennb,'r','r')
          call rdedx(q(kfao),lennb,kvb7_fcao,num8)
          call mult11(q(kfao),q(kfmo),q(kscr),nn,nbasis,q(kvec),
     &                q(kdag))
          if (iprint.gt.1500) then
            write(iwr,*) 'fock matrix on mo basis (including core)'
            call tripri(q(kfmo),nn)
          end if
c...  Store fock matrix on mo basis on ed7
          kvb7_fcmo = kscra7vb('kvb7_fcmo',lennn,'r','w')
          call wrt3(q(kfmo),lennn,kvb7_fcmo,num8)
        endif
        timing = cpulft(1) - timing
      endif
c
c...  set ig(1 to 4,0) to 0, and ig(5,0) to 1, as the number of alpha
c...  blocks might be zero (ialfa=zero, see symblo). ig is used starting
c...  from the sixth position ! (see call hamilt)
c
      ig(1) = 0
      ig(2) = 0
      ig(3) = 0
      ig(4) = 0
      ig(5) = 1
c
c...  calculate the matrix-elements
c
      nstruc = nbrill
      kgrh = 1
      kgrs = kgrh + nbrill
      kdiagh = kgrs + nbrill
      kdiags = kdiagh + nbrill
      nused = kdiags + nbrill
      kir = nused
      kic = kir + (nbrill-1)/nipw() + 1 
      nused = kic + (nbrill-1)/nipw() + 1
      nword = lword - nused
      if (nused.gt.lword) call corfait(nused,lword,'before call bamilt')
      timeh = cpulft(1)
      call bamilt(pacdet,idet,jdet,detcomb,dettot,icp,jcp,
     &            idetps,coeff(ncoeff+1),ig(6),igroup,ndetps,
     &            hamil,overl,q(kic),q(kir),scr1,scr2,scr3,s,g,
     &            supers,superh,superg(2),ipos,weight,
     &            q(kgrh),q(kgrs),q(kdiagh),q(kdiags),
     &            nelec,ngroup,nalfa,iortho,northo,q(nused),nword)
      timeh = cpulft(1) - timeh
      if (iprinv.gt.10) write(iwr,66) timeh
66    format(/' ','the construction of the matrix representation costs '
     &       ,f8.3,' seconds')
      return
      end

c***********************************************************************
      subroutine breorb(idet,ndets,nelec,nalfa,iorth,coeff)
c
      implicit REAL  (a-h,o-z) , integer   (i-n)
c
c
c...  "reorb" for brillouin functions : is slightly simpler because
c...  in the brillouin function determinants only occur once => no
c...  scratch.
c...  in this routine the determinants in idet are reordered so that
c...  the orbitals appear in mutually orthogonal groups, for the benifit
c...  of the matrix element evaluation. this is done for alpha and beta
c...  spin separately, as spin introduces an extra orthogonality on top
c...  of the spatial one. in the i-th position iorth must contain a
c...  number between 1 and the number of orthogonality groups, that is
c...  characteristic for the i-th orbital its group. parity changes are
c...  accounted for via the coefficients of the determinants.
c
      dimension idet(nelec,ndets),iorth(*),coeff(ndets)
      do m=1,ndets
         i     =  1
         irrep = 1
10       do k=i,nalfa
            if (irrep.eq.1) then
               if (iorth(idet(k,m)).eq.0) call caserr('boem joop')
            end if
           if (irrep.eq.iorth(idet(k,m))) then
               if (i.ne.k) then
                  ipar = -ipar
                  iii        =  idet(k,m)
                  idet(k,m)  =  idet(i,m)
                  idet(i,m)  =  iii
                  coeff(m) = -coeff(m)
               end if
               i = i + 1
            end if
         end do
         if (i.le.nalfa) then
            irrep = irrep + 1
            goto 10
         end if
         irrep  = 1
30       do k=i,nelec
            if (irrep.eq.1) then
               if (iorth(idet(k,m)).eq.0) call caserr('boem joop beta')
            end if
            if (irrep.eq.iorth(idet(k,m))) then
               if (i.ne.k) then
                  iii       =  idet(k,m)
                  idet(k,m) =  idet(i,m)
                  idet(i,m) =  iii
                  coeff(m) = -coeff(m)
               end if
               i = i + 1
            end if
         end do
         if (i.le.nelec) then
            irrep = irrep + 1
            goto 30
         end if
      end do
c
      return
      end

c***********************************************************************
      subroutine brigen(nelec,nalfa,igroup,ngroup,pacdet,
     &                  coeff,bcoeff,ciref,ndetps,idetps,idet,
     &                  ibdet,iortho,northo,q,lword,superh,supers,
     &                  superg)
c
      implicit REAL (a-h,o-z), integer (i-n)
c
INCLUDE(common/vbcri)
INCLUDE(common/scftvb)
INCLUDE(../m4/common/sizes)
INCLUDE(../m4/common/iofile)
INCLUDE(../m4/common/files)
INCLUDE(common/ffile)
INCLUDE(common/turtleparam)
INCLUDE(common/c8_16vb)
INCLUDE(common/twice)
INCLUDE(common/splice)
INCLUDE(common/tractlt)
INCLUDE(common/brill)
INCLUDE(common/vbpert)
      common/blkin/gg(340),mij(340),mword,mdumm
_IF(linux)
      external fget
_ENDIF
      logical notnull
      dimension igroup(5,*),pacdet(*),coeff(*),ciref(*),ibdet(nelec,*),
     &          ndetps(*),idetps(*),idet(nelec,*),bcoeff(*),iortho(*),
     &          q(*),superh(*),superg(*),supers(*)
c
      nvarp = 1
      ndet = igroup(1,1)
      nwp = (ndet*nelec-1) / (64/n8_16) + 1
_IF(atmol)
      rewind 33
c...  read "new" coefficients
      read(33) ndum,(bcoeff,ij=kbold,ndet)
_ELSE
      kvb7_bcn = kscra7vb('kvb7_bcn',igroup(1,1),'r','n')
      if (kvb7_bcn.lt.0) then
         call caserr('Psi0 coefs needed in brigen but not found on ed7')
      else
         kvb7_bcn = kscra7vb('kvb7_bcn',igroup(1,1),'r','r')
         call rdedx(bcoeff,igroup(1,1),kvb7_bcn,num8)
      endif
_ENDIF
_IF(atmol)
         rewind 25
         read(25) nelec,nalfa,nconf,ndets,nstruc,maxdet,maxstr,
     &         nwpack,ncoeff,imax
_ELSE
         call wr15vb(nelec,nalfa,nconf,ndets,nstruc,maxdet,maxstr,
     &               nwpack,ncoeff,imax,ndum11,ndum12,'read')
_ENDIF
      ngroup      = 1
      it = ndet
      kkq = 1
c
      do 100 iq=1,nequi
c
c...  loop over groups of occupieds to be changed
c
         nrejec = 0
         jend = nex(kkq)+1
c
         do 90 j=2,jend
            jj = j - nrejec
            ngroup = ngroup + 1
            lngroup= 0
            nb = 0
            is = it
            do 85 jq=1,iequi(iq)
               ibaal = kkq+jq-1
               ifr = iex(1,ibaal)
               ito = iex(jj,ibaal)
               ione = 1
               signn = + 1
79             continue
               do 80 k=1,ndet
                  call bcreat(idet(1,k),nelec,nalfa,ibdet(1,nb+1),nbb,
     &                                                        ifr,ito)
c...  nbb is # new brillouin dets.,resulting out of the
c...  reference det. idet(1,k), from the excitation ifr=>ito
                  l = nb
77                if (l.lt.nb+nbb) then
c...  new dets. have to be checked
                     l = l + 1
                     do 65 m=1,nb
c...  see if det. l occured before
                        isignq = isame(ibdet(1,m),ibdet(1,l),
     &                             ibdet(1,nb+nbb+1),nelec,nalfa)
                        if (isignq.ne.0) then
c...  det. l occurred before, put it onto the
c...  older one (that is: add the coefficients)
                           bcoeff(is+m) = bcoeff(is+m)
     &                             + isignq * bcoeff(k) * iflip(ifr,ito)
                           l = l - 1
c...  now forget about det. l
                           nbb = nbb - 1
                           if (l.lt.nb+nbb) then
c...  copie the next one onto the one just
c...  thrown away
                              do 63 n=1,nelec
c63 bug (cb 16-08-94)         ibdet(n,l)   = ibdet(n,l+1)
 63                           ibdet(n,l+1) = ibdet(n,l+2)
                           end if
                           go to 77
                        end if
65                   continue
                     it = it + 1
                     bcoeff(it) = bcoeff(k) * iflip(ifr,ito) * signn
                     go to 77
                  end if
                  nb = nb + nbb
80             continue
*
               if (unitary.and.ione.eq.1) then
                  ione = 2
                  ifr = ito
                  ito = iex(1,ibaal)
                  signn = - 1
                  go to 79
               end if
*
85          continue
            notnull = .false.
            do 87 k=is+1,it
               if (dabs(bcoeff(k)).gt.cribri) notnull = .true.
87          continue
            if (.not.notnull) then
c
c...  symmetry or our intervention
c...  causes annihilation of states (see above)
c
               nb = 0
               it = is
c...  print includes core in numbering
               write(iwr,602)
     &         (iex(1,m),iex(jj,m)+ncore,m=kkq,kkq+iequi(iq)-1)
602            format(' annihilation of ',(t17,5(i4,' =>',i4)))
            end if
            if (nb.eq.0) then
c...  no new dets => remove excitation(s)
               if (ovbper.and.igroup(4,ngroup).eq.1) then
                  nvarp = nvarp - 1
               end if
               ncopy = nex(kkq)-jj+1
               do 89 k=kkq,kkq+iequi(iq)-1
                  do l=1,ncopy
                     iex(jj+l-1,k) = iex(jj+l,k) 
                     ieqsig(jj+l-1,k) = ieqsig(jj+l,k) 
                  end do
                  nex(k) = nex(k) - 1
89             continue
               nrejec = nrejec + 1
               ngroup = ngroup - 1
               goto 90
            end if
c
c...  reorder determinants so as to group the orbitals into
c...  orthogonality groups. this is vital for matre3 (otherwise
c...  it produces nonsense). match tries to find per det its
c...  complement (in terms of alpha-beta) see call to mathab
c
            call breorb(ibdet(1,1),nb,nelec,nalfa,iortho,bcoeff(is+1))
            call match(ibdet,ibdet(1,nb+1),nb,nelec,nalfa,
     &                        ndetps(ngroup),idetps(is+1))
            call pack(pacdet(nwp+1),n8_16,ibdet(1,1),nb*nelec)
            nwp = nwp + (nb*nelec-1) / (64/n8_16) + 1
            igroup(1,ngroup) = nb
            igroup(2,ngroup) = iq
            igroup(3,ngroup) = nb
90       continue
         nbrill=ngroup
100   kkq = kkq + iequi(iq)
      return
      end

c***********************************************************************
      integer function isame(idet1,idet2,icp,nelec,nalfa)
c
c...  check if idet1 and idet2 are build from the same integers
c     value returned will be parity if so
c
      implicit REAL (a-h,o-z)
      dimension idet1(nelec),idet2(nelec),icp(nelec)
c
c...  check a simple criterion
c
      isame = 0
      isum1 = idet1(1)
      isum2 = idet2(1)
      do 1 i=2,nalfa
         isum1 = isum1 + idet1(i)
         isum2 = isum2 + idet2(i)
1     continue
      if (isum1.ne.isum2) return
      do 2 i=nalfa+1,nelec
         isum1 = isum1 + idet1(i)
         isum2 = isum2 + idet2(i)
2     continue
      if (isum1.ne.isum2) return
      iperm = 1
      do 10 i=1,nelec
10    icp(i) = idet2(i)
      do 100 ir=1,nalfa-1
      if (idet1(ir).ne.icp(ir)) then
c...  look for idet1(ir)
         do 20 ic=ir+1,nalfa
20       if (icp(ic).eq.idet1(ir)) go to 30
         isame = 0
c...  not found, so not the same
         return
c...  interchange ic and ir
30       iperm = -iperm
         isave = icp(ir)
         icp(ir) = icp(ic)
         icp(ic) = isave
      end if
100   continue
      do 200 ir=nalfa+1,nelec-1
      if (idet1(ir).ne.icp(ir)) then
c...  look for idet1
         do 120 ic=ir+1,nelec
120      if (icp(ic).eq.idet1(ir)) go to 130
         isame = 0
c...  not found, so not the same
         return
c...  interchange ic and ir
130      iperm = -iperm
         isave = icp(ir)
         icp(ir) = icp(ic)
         icp(ic) = isave
      end if
200   continue
c
      isame = iperm
c
      return
      end

c***********************************************************************
      integer function isameb(idet1,idet2,icp,nelec,nalfa)
      implicit REAL  (a-h,o-z) ,  integer   (i-n)
c
c...  almost the same as isame, this one compairs alpha orbitals
c...  of idet1 with beta orbitals of idet2 (for "match")
c...  check if idet1 and idet2 are build from the same integers
c...  value returned will be parity if so
c...  note that nalfa must equal nbeta
c
      dimension idet1(nelec),idet2(nelec),icp(nelec)
c
c...  check a simple criterion
c
      if(nelec-nalfa.ne.nalfa) then
         isameb = 0
         return
      end if
      isameb = 0
      isum1 = idet1(1)
      isum2 = idet2(nalfa+1)
      do 1 i=2,nalfa
         isum1 = isum1 + idet1(i)
         isum2 = isum2 + idet2(nalfa+i)
1     continue
      if (isum1.ne.isum2) return
      do 2 i=nalfa+1,nelec
         isum1 = isum1 + idet1(i)
         isum2 = isum2 + idet2(i-nalfa)
2     continue
      if (isum1.ne.isum2) return
      iperm = 1
      do 10 i=1,nalfa
10    icp(i) = idet2(i+nalfa)
      do 12 i=nalfa+1,nelec
12    icp(i) = idet2(i-nalfa)
      do 100 ir=1,nalfa-1
      if (idet1(ir).ne.icp(ir)) then
c...  look for idet1(ir)
         do 20 ic=ir+1,nalfa
20       if (icp(ic).eq.idet1(ir)) go to 30
         isameb = 0
c...  not found, so not the same
         return
c...  interchange ic and ir
30       iperm = -iperm
         isave = icp(ir)
         icp(ir) = icp(ic)
         icp(ic) = isave
      end if
100   continue
      do 200 ir=nalfa+1,nelec-1
      if (idet1(ir).ne.icp(ir)) then
c...  look for idet1
         do 120 ic=ir+1,nelec
120      if (icp(ic).eq.idet1(ir)) go to 130
         isameb = 0
c...   not found, so not the same
         return
c...  interchange ic and ir
130      iperm = -iperm
         isave = icp(ir)
         icp(ir) = icp(ic)
         icp(ic) = isave
      end if
200   continue
c
      isameb = iperm
c
      return
      end

c***********************************************************************
      subroutine bsmat(supers,superh,superg,ciref,nbasis,nsa,nvirt)
      implicit REAL (a-h,o-z), integer (i-n)
INCLUDE(../m4/common/sizes)
INCLUDE(../m4/common/vcore)
INCLUDE(common/turtleparam)
INCLUDE(common/c8_16vb)
c
c...  dynamical core allocation and brillouin control routine
c
c... for dynamic memory allocation debugging purpose, contains IGMEM_ vars.
INCLUDE(../m4/common/gmemdata)
INCLUDE(../m4/common/gmempara) 
c  
      dimension supers(*),superh(*),superg(*),ciref(*)
c
      common /icaselist/ icaselist(35)
      common /array/ narray,iarray(100)
      common /arnam/ anames(100)
INCLUDE(common/brill)
INCLUDE(common/hsinfo)
c
      character*8 anames
c
_IF(atmol)
      rewind 25
      read(25) nelec,nalfa,nconf,ndets,nstruc,maxdet,maxstr,
     &         nwpack,ncoeff,imax
_ELSE
      call wr15vb(nelec,nalfa,nconf,ndets,nstruc,maxdet,maxstr,
     &            nwpack,ncoeff,imax,ndum11,ndum12,'read')
_ENDIF
c      do ii=1,35
c        icaselist(ii) = 0
c      end do
       call putintdev(superg)
      npword = nwpack
      nextot = 0
      maxequ = 0
      do 5 i=1,nequi
5     maxequ = max(iequi(i),maxequ)
      maxdet = max(ndets,2*ncoeff*maxequ)
      do 10 i=1,nsa
         do 15 j=1,nex(i)
            npword = npword + (2*ndets*nelec-1) / (64/n8_16) + 1
15       continue
         nextot = nextot + nex(i)
10    continue
      nbdets = 2 * ncoeff * nextot + ncoeff
      narray = 0
c
c...  pacdet will contain the packed slater determinants. these are
c...  packed per "brillouin-structure" (bs)
c
      kpacde = igmem_alloc_inf(npword,'vbscf.m','bsmat','kpacde',
     &                         IGMEM_DEBUG)
c
c...  idet and jdet contain the unpacked slater determinants per bs
c...  in "brigen" : idet will contain all the determinants that make up
c...  the reference function, jdet will contain the brillouin
c...  determinants that arise because of one excitation  + 2 scratch det
c
      kidet  = igmem_alloc_inf(maxdet*nelec,'vbscf.m','bsmat','kidet',
     &                         IGMEM_DEBUG)
      kjdet  = igmem_alloc_inf((maxdet+2)*nelec,'vbscf.m','bsmat',
     &                         'kjdet',IGMEM_DEBUG)
c
c...  detcomb is scratch for the calculation of the brill. matrix-
c...  elements, dettot is scratch for the corresponding overlaps
c
      kdetco = igmem_alloc_inf(maxdet,'vbscf.m','bsmat','kdetco',
     &                         IGMEM_DEBUG)
      kdetto = igmem_alloc_inf(maxdet,'vbscf.m','bsmat','kdetto',
     &                         IGMEM_DEBUG)
c
c...  the matrix element routine (matre3) is given the determinants via
c...  icp and jcp (mind parity changes !).
c
      kicp   = igmem_alloc_inf(nelec,'vbscf.m','bsmat','kicp',
     &                         IGMEM_DEBUG)
      kjcp   = igmem_alloc_inf(nelec,'vbscf.m','bsmat','kjcp',
     &                         IGMEM_DEBUG)
c
c...  ndetps contains the number of determinants per structure
c
      kndetp = igmem_alloc_inf(nstruc+max((nextot+1),nstruc),'vbscf.m',
     &                         'bsmat','kndetp',IGMEM_DEBUG)
c
c...  idetps contains the determinant numbers that define the structures
c...  coeff contains the corresponding coefficients + ncoeff /scratch
c
      kidetp = igmem_alloc_inf(nbdets,'vbscf.m','bsmat','kidetp',
     &                         IGMEM_DEBUG)
      kncoef = igmem_alloc_inf(nbdets+ncoeff,'vbscf.m','bsmat','kncoef',
     &                         IGMEM_DEBUG)
c
c...  ig will contain the blocking information per matrix element
c...  (5 numbers per block). but at first (in symblo) it contains
c...  5 numbers per orthogonality group => at most norb * 10 numbers !
c...  (spin-orthogonality on top of spatial orthogonality => 2 * 5)
c...  ig(1 to 4,0) is put to zero, ig(5,0)=1, this must be done to deal
c...  with the case no alpha block occur in the overlap matrix (ialfa=0)
c...   => at most norb * 10 + 5 numbers

      kig    = igmem_alloc_inf(10*(nsa+nvirt)+5,'vbscf.m','bsmat','kig',
     &                         IGMEM_DEBUG)
c
c...  igroup contains the number of determinants per bs, the number
c...  of structures per bs (=1) and the starting position of idetps per
c...  group.
c
      kigrou = igmem_alloc_inf(max(nconf,(nextot+1))*5,'vbscf.m',
     &                         'bsmat','kigrou',IGMEM_DEBUG)
c
c...  hamil and overl contain the matrix representation of the
c...  hamiltonian on a structure basis
c
      khamil = igmem_alloc_inf(nextot+1,'vbscf.m','bsmat','khamil',
     &                         IGMEM_DEBUG)
      koverl = igmem_alloc_inf(nextot+1,'vbscf.m','bsmat','koverl',
     &                         IGMEM_DEBUG)
c
c...  ipos is the position array that is used to gather the integrals
c...  per matrix element. it is also used to reorder the matrix
c...  elements before they are written to disc (in writh)
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
c...
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
      kipos  = igmem_alloc_inf(max(nnnnn,nextot+1),'vbscf.m','bsmat',
     &                         'kipos',IGMEM_DEBUG)
      kweigh = igmem_alloc_inf(nwwww,'vbscf.m','bsmat','kweigh',
     &                         IGMEM_DEBUG)
c
c...  g contains the integrals that make up the matrix element
c
      kg     = igmem_alloc_inf(nnnnn,'vbscf.m','bsmat','kg',IGMEM_DEBUG)
c
c...  iortho contains the orthogonality number per orbital
c
      kiorth = igmem_alloc_inf(nsa+nvirt,'vbscf.m','bsmat','kiorth',
     &                         IGMEM_DEBUG)
c
c...  s contains the (biorthogonalised) overlap-matrix per matrix-
c...  element
c
      ks     = igmem_alloc_inf(nalfa**2+nbeta**2,'vbscf.m','bsmat','ks',
     &                         IGMEM_DEBUG)
c
c...  scr1 and scr2 are used in hamilt during the transformation of
c...  hamil and overl. scr1/scr2/scr3 are used in matre3 using
c...  always less than (nelec**2) words. scr1 is used in hamilt to
c...  gather (reorder) matrix elements before they are written to disc.
c
      nnnnn  = max( nelec**2 , maxdet * 1 )
      mmmmm  = max( nnnnn    , nextot + 1 )
      kscr1  = igmem_alloc_inf(mmmmm,'vbscf.m','bsmat','kscr1',
     &                         IGMEM_DEBUG)
      kscr2  = igmem_alloc_inf(nnnnn,'vbscf.m','bsmat','kscr2',
     &                         IGMEM_DEBUG)
      kscr3  = igmem_alloc_inf(nelec*nelec,'vbscf.m','bsmat','kscr3',
     &                         IGMEM_DEBUG)
      kmaxmem = igmem_max_memory()
      kmemscr = kmaxmem/10
      kscr   = igmem_alloc_inf(kmemscr,'vbscf.m','bsmat','kscr',
     &                         IGMEM_DEBUG)
      call bmain(Q(kpacde),Q(kidet),Q(kjdet),
     &           Q(kdetco),Q(kdetto),
     &           Q(kicp),Q(kjcp),Q(kidetp),
     +           Q(kncoef),Q(kig),
     &           Q(kigrou),Q(khamil),Q(koverl),
     &           supers,superh,superg,
     +           Q(kipos),Q(kweigh),
     &           Q(kg),Q(kiorth),Q(ks),
     +           Q(kscr1),Q(kscr2),
     &           Q(kscr3),Q(kndetp),
     &           ciref,Q(kscr),kmemscr,nvirt)

c
c      call bmain(qq(kpacde),qq(kidet),qq(kjdet),qq(kdetco),qq(kdetto),
c     &           qq(kicp),qq(kjcp),qq(kidetp),qq(kncoef),qq(kig),
c     &           qq(kigrou),qq(khamil),qq(koverl),
c     &           supers,superh,superg,qq(kipos),qq(kweigh),
c     &           qq(kg),qq(kiorth),qq(ks),qq(kscr1),qq(kscr2),
c     &           qq(kscr3),qq(kndetp),
c     &           ciref,qq(kscr),kmemscr,nvirt)
c
      call gmem_free_inf(kscr,'vbscf.m','bsmat','kscr')
      call gmem_free_inf(kscr3,'vbscf.m','bsmat','kscr3')
      call gmem_free_inf(kscr2,'vbscf.m','bsmat','kscr2')
      call gmem_free_inf(kscr1,'vbscf.m','bsmat','kscr1')
      call gmem_free_inf(ks,'vbscf.m','bsmat','ks')
      call gmem_free_inf(kiorth,'vbscf.m','bsmat','kiorth')
      call gmem_free_inf(kg,'vbscf.m','bsmat','kg')
      call gmem_free_inf(kweigh,'vbscf.m','bsmat','kweigh')
      call gmem_free_inf(kipos,'vbscf.m','bsmat','kipos')
      call gmem_free_inf(koverl,'vbscf.m','bsmat','koverl')
      call gmem_free_inf(khamil,'vbscf.m','bsmat','khamil')
      call gmem_free_inf(kigrou,'vbscf.m','bsmat','kigrou')
      call gmem_free_inf(kig,'vbscf.m','bsmat','kig')
      call gmem_free_inf(kncoef,'vbscf.m','bsmat','kncoef')
      call gmem_free_inf(kidetp,'vbscf.m','bsmat','kidetp')
      call gmem_free_inf(kndetp,'vbscf.m','bsmat','kndetp')
      call gmem_free_inf(kjcp,'vbscf.m','bsmat','kjcp')
      call gmem_free_inf(kicp,'vbscf.m','bsmat','kicp')
      call gmem_free_inf(kdetto,'vbscf.m','bsmat','kdetto')
      call gmem_free_inf(kdetco,'vbscf.m','bsmat','kdetco')
      call gmem_free_inf(kjdet,'vbscf.m','bsmat','kjdet')
      call gmem_free_inf(kidet,'vbscf.m','bsmat','kidet')
      call gmem_free_inf(kpacde,'vbscf.m','bsmat','kpacde')

      call delintdev(superg)
      return
      end

c***********************************************************************
      subroutine determ(s,ndim,det,irank)
      implicit REAL  (a-h,o-z) , integer   (i-n)
c.....
c.....this routine determines the determinant of s. s is spoiled
c.....
INCLUDE(common/vbcri)
      dimension s(ndim,ndim)
      irank = 0
      det   = 1.0d0
      if (ndim.eq.0) return
      ipar  = 1
      call piv(s,ndim,ndim,1,pivo,ipar)
      if (dabs(pivo).lt.cridep) then
         irank = 0
         det   = 0.0d0
         return
      end if
      irank = 1
      do 40 i = 1,ndim-1
         det = det * s(i,i)
         do 30 j = i+1,ndim
c.....
c.....form i-th row of r
c.....
            s(i,j) = -s(i,j)/s(i,i)
c.....
c.....adapt s-matrix to r
c.....
            do 10 k = i+1,ndim
               s(k,j) = s(k,j) + s(k,i)*s(i,j)
10          continue
c.....
c.....adapt r-matrix to itself
c.....
            do 20 l = 1,i-1
               s(l,j) = s(l,j) + s(l,i)*s(i,j)
20          continue
30       continue
         call piv(s,ndim,ndim,i+1,pivo,ipar)
         if (dabs(pivo).lt.cridep) then
            det = 0.0d0
            return
         end if
         irank = irank + 1
40    continue
      det = det * s(ndim,ndim) * ipar
      return
      end

c***********************************************************************
      subroutine vbdiis(ro,vo,vn,nbasis,nact,cr,lword,maxdim,mindim)
c
      implicit REAL (a-h,o-z), integer(i-n)
c
c...  subroutine to perform DIIS extrapolation for VBSCF
c...
c...  DIIS tries to find linear combination of nit vectors
c...  that minimises the corresponding error-vectors in a least
c...  square-norm . (note nact mo's is considered 1 vector)
c...
c...  see : P.Pulay  Chem.Phys. Lett. 73, 393 (1980)
c...        P.Pulay  J.Comput.Chem. 3, 556 (1982)
c...        T.P. Hamilton, P.Pulay J.Chem.Phys. 84, 5728 (1986)
c...        J.H. van Lenthe, J. Verbeek, P.Pulay  Mol. Phys. 73(5), 1159-1170. (1991)
c
c
c...   this version uses as error-vectors the differences between
c...   new and old-vectors in orbopt and as vectors the new-vectors
c...
c...   a future version might use the orbital-rotation gradients instead
c...   and so avoid building the whole BI-matrix
c...   Common basis :  symmetrically orthogonalised AO's
c
c...   called from orbopt
c
c...   MANUAL : (as implemented in scfin)
c...   diis iwdiis ntdiis madiis midiis
c...     where :
c...     iwdiis/ntdiis  : if "del" (log(delta psi)) < iwdiis   and
c...                      # iterations >= ntdiis : start DIIS
c...     madiis : Max. # DIIS iterands (default 20)
c...     midiis : Min. # DIIS iterands in interpolation (default 3)
c...   madiis and midiis arethe formal parameters  maxdim and mindim
c...   just giving diis in the scf-input (so del < 0 and no restriction
c...   on # iterations) seems fine
c
c...  File-usage :
c...    unit 1 (iun1) : old vectors
c...    unit 2 (iun2) : old error-vectors
c...    unit 3 (iun3) : transformation-matrix/DIIS-matrix
c...   They are in GAMESS verions parts of ed7
c...  It keeps track of the times it's been called
c...  *** it is assumed that no sudden phase changes occur ***
c...      ORBOPT has been changed to (hopefully) ensure this by
c...      making sure that BI(1) > 0
c...      this is checked in DIISPH, where also normalisations occur
c
c...  vo(nbasis,nact) : old orbitals
c...  vn(nbasis,nact) : new orbitals
c...  cr             : scratch-space (lword)
c...  maxdim : maximum dimension of diis space
c...
c...  Joop van Lenthe , Fayetteville, Arkansas 1989
c...  first unreliable etc version ... seems to work though
c
c///   diirw is used because compiler has trouble with implied do
c
      dimension vo(nbasis,nact),vn(nbasis,nact),cr(lword),ro(*)
INCLUDE(common/scftvb)
INCLUDE(common/vbdiist)
INCLUDE(common/vbcri)
INCLUDE(../m4/common/iofile)
INCLUDE(common/scra7vb)
c
      parameter (iun1=1,iun2=2,iun3=3)
c...   1,2,3, are kept to indicate vb1,vb2,vb3 also on ed7
c
      save nit
c
c...  determine current diis dimension
c...  if larger than current # iterations first few vectors/gradients
c...  must be skipped
c
      if (maxdim.lt.3) call vberr('no diis  < 3')
      nit = nit + 1
      ndiim = min(nit,maxdim)
      nskip = nit - ndiim
c
c     nskip is vectors to be skipped in atmol version
c     in gamess version we use a mod(i,maxdim=maxdis)
c
c...  clear diagonal (we want do make gradient small)
c...  keep diagonal / we may want to put it back
c
      ktr = 1
c
      if (nit.eq.1) then
_IFN(atmol)
c
c...  layout of VB1,VB2,VB3 on ed7 
c...    unit 1 (VB1) : old vectors   - maxdim * nbasis*nact
c...    unit 2 (VB2) : old error-vectors - maxdim * nbasis*nact
c...    unit 3 (VB3) : transformation-matrix/DIIS-matrix : 2 * nbasis*nbasis, maxdim*(maxdim+1)/2
c      This is arranged HERE
c
         kvb1 = kscra7vb('kvb1',0,'r','n')
         if (kvb1.gt.0.and.maxdis.ne.maxdim) 
     1      call caserr('diis consistency error')
         maxdis = maxdim
         nvb12 = lensec(nbasis*nact)
         kvb1 = kscra7vb('kvb1',maxdim*nvb12,'b','w')
         kvb2 = kscra7vb('kvb2',maxdim*nvb12,'b','w')
         nvb3 = lensec(nbasis*nbasis)
         nvb3a = lensec(maxdim*(maxdim+1)/2)
         kvb3 = kscra7vb('kvb3',2*nvb3+nvb3a,'b','w')
c        kbl7vb = kvb3 + 2*nvb3 + lensec(maxdim*(maxdim+1)/2)
_ENDIF
c
c...     if first time  // determine the transformation to
c...     orthogonal ao's
c
         ks = ktr + nbasis*nbasis
         kcr = ks + nbasis*(nbasis+1)/2
         if (kcr+nbasis.gt.lword) call caserr('core fail1 diis')
c
c...     get s-matrix
c
         call get1e(cr(ks),dummy,'s',cr(ks))
c
         call sminh(cr(ks),cr(ktr),nbasis,cr(kcr))
c
c...     in cr(ktr) transformation from AO's to orthogonal AO's
c...     open vector-file and dump transformation to it
c
c         maxrec = max(nbasis*nbasis,nbasis*nact)*8
c...CB   recl is only valid for direct-access files
c        open(iun1,form='unformatted',status='scratch',
c    *        access='sequential',recl=200000)
c        open(iun2,form='unformatted',status='scratch',
c    *        access='sequential',recl=200000)
c        open(iun3,form='unformatted',status='scratch',
c    *        access='sequential',recl=200000)
c        open(iun1,form='unformatted',status='scratch',
c    *        access='sequential')
c        open(iun2,form='unformatted',status='scratch',
c    *        access='sequential')
c        open(iun3,form='unformatted',status='scratch',
c    *        access='sequential')
c
_IF(atmol)
         call fopen(iun1,'VB1','priv','unformatted')
         call fopen(iun2,'VB2','priv','unformatted')
         call fopen(iun3,'VB3','priv','unformatted')
         call diirw(cr(ktr),nbasis*nbasis,iun3,'write')
_ELSE
         kvb3 = kscra7vb('kvb3',nbasis*nbasis,'r','n')
         call wrt3(cr(ktr),nbasis*nbasis,kvb3,num8)
_ENDIF
c
c...     invert the matrix and dump
c
         call osinv_vb(cr(ktr),nbasis,determ,crilow,cr(ks),cr(kcr))
c
c...     now in cr(ktr) transformation from orthog AO's to AO's
c
_IF(atmol)
         call diirw(cr(ktr),nbasis*nbasis,iun3,'write')
_ELSE
         kvb3 = kscra7vb('kvb3',nbasis*nbasis,'r','n')
         call wrt3s(cr(ktr),nbasis*nbasis,num8)
_ENDIF
c
      else
c
c...   read transformation-matrix back  (need inverse first)
c
_IF(atmol)
         rewind iun3
         read(iun3)
         call diirw(cr(ktr),nbasis*nbasis,iun3,'read')
_ELSE
         kvb3 = kscra7vb('kvb3',nbasis*nbasis,'r','n')
         call rdedx(cr(ktr),nbasis*nbasis,
     1              kvb3+lensec(nbasis*nbasis),num8)
_ENDIF
c
      end if
c
      kvo = ktr + nbasis*max(nbasis,nact)
      kvn = kvo + nbasis*nact
      if (kvn+nbasis*nact.gt.lword) call caserr('core fail2 diis')
c
c...  transform old and new vectors to orthonormal basis
c...  using orthog AO => AO transformation matrix
c
      call mxmd(cr(ktr),1,nbasis,vo,1,nbasis,cr(kvo),1,nbasis,
     *          nbasis,nbasis,nact)
c
      call mxmd(cr(ktr),1,nbasis,vn,1,nbasis,cr(kvn),1,nbasis,
     *          nbasis,nbasis,nact)
c      kvnn = kvn + nbasis*nact
c      call fmove(ro,cr(kvnn),nbasis*nact)
c      call mxmd(cr(ktr),1,nbasis,cr(kvnn),1,nbasis,ro,1,nbasis,
c     &          nbasis,nbasis,nact)
c
c...  check phase of vectors against the previous ones (ktr=scratch)
c...  and write current vectors
c...  also check new agains old vectors
c...  and normalise
c
      if (nit.gt.1) then
_IF(atmol)
         rewind iun1
         do 10 i=1,nit-2
10       read(iun1)
         call diirw(cr(ktr),nbasis*nact,iun1,'read')
_ELSE
         kkd = mod(nit-2,maxdis)
         kvb1 = kscra7vb('kvb1',0,'r','n')
         call rdedx(cr(ktr),nbasis*nact,kvb1+kkd*nvb12,num8)
_ENDIF
         call diisph(cr(ktr),cr(kvo),nbasis,nact,'previous vs old')
      end if
      call diisph(cr(kvo),cr(kvn),nbasis,nact,'old vs new')
_IF(atmol)
      call diirw(cr(kvn),nbasis*nact,iun1,'write')
_ELSE
      kkd = mod(nit-1,maxdis)
      kvb1 = kscra7vb('kvb1',0,'r','n')
      call wrt3(cr(kvn),nbasis*nact,kvb1+kkd*nvb12,num8)
_ENDIF
c
c...  compute error-vector vn-vo  in cr(ktr)
c
      call subvec(cr(ktr),cr(kvn),cr(kvo),nbasis*nact)
c      call fmove(ro,cr(ktr),nbasis*nact)
c
c...  set up diis matrix <error!error>  (triangle)
c...  diism_vb reads old matrix and dumps the new one
c...  (assuming only 2 record before it on iun3)
c...  diism_vb dumps the last error-vector as well
c
      kdm = ktr + nbasis*max(nact,nbasis)
      kcr = kdm + ndiim*(ndiim+1)/2
      kend = kcr + nbasis*nact
      if (kend.gt.lword) call vberr('core-failure 2 in vbdiis')
c
      call diism_vb(cr(ktr),cr(kcr),nbasis*nact,cr(kdm),ndiim,
     *           iun2,nskip,iun3)
c
c...  if insufficient vectors gathered return
c
      if (ndiim.lt.mindim) then
         call fmove(vn,vo,nbasis*nact)
         ovbdiis = .false.
         return
      end if
c
c...  solve diis equations
c
      kres = kdm
      kdmn = kres + ndiim*(ndiim+1)/2
      kscr1 = kdmn + ndiim*ndiim
      kscr2 = kscr1 + ndiim
      kend = kscr2 + ndiim
      if (kend.gt.lword) call vberr('core-failure 3 in vbdiis')
c
      call diiseq(cr(kdmn),cr(kscr1),cr(kscr2),ndiim,deter,cr(kdm))
c
c...  print solution for a while  (**debug**)
c
      rmax = 0.0d0
      maxdiis = 0
      abssum = 0.0d0
      do 20 i=1,ndiim
         if (dabs(cr(kres+i-1)).gt.rmax) then
            rmax = dabs(cr(kres+i-1))
            maxdiis = i
         end if
         abssum = abssum + dabs(cr(kres+i-1))
20    continue
      percmax = 100.d0 * rmax / abssum
      perccur = 100.d0 * dabs(cr(kres+ndiim-1)) / abssum
      nitdiis = ndiim
      if (iprins.gt.1) write(iwr,988) deter,(cr(kres-1+kk),kk=1,ndiim)
988    format(5x,'DIIS det ',e11.4,' coeff ...',/,(5x,8F13.6))
c
c...  determine extrapolated vectors
c
      kcr = ktr
      kv = kres + ndiim
      kend = kv + nbasis*nact
      if (kend.gt.lword) call vberr('core-failure 4 in vbdiis')
      call diisrs(cr(kv),cr(kcr),nbasis*nact,cr(kres),ndiim,iun1,nskip)
c
c...  that's it
c..   now transform vectors to ao-basis
c...  use AO to orthog AO now
c
_IF(atmol)
      rewind iun3
      call diirw(cr(ktr),nbasis*nbasis,iun3,'read')
_ELSE
      kvb3 = kscra7vb('kvb3',0,'r','n')
      call rdedx(cr(ktr),nbasis*nbasis,kvb3,num8)
_ENDIF
      call mxmd(cr(ktr),1,nbasis,cr(kv),1,nbasis,vo,1,nbasis,
     *          nbasis,nbasis,nact)
c
      return
c
      entry init_vbdiis
c
      nit = 0
      nitdiis = 0
c     
      return
      end

c***********************************************************************
_IF(atmol)
      subroutine diirw(q,n,iun,action)
c
c...   avoid implied do-loop problems
c
      implicit REAL (a-h,o-z), integer (i-n)
      character*(*) action
      dimension q(n)
c
      if (action.eq.'read') then
         read(iun) q
      else if (action.eq.'write') then
         write(iun) q
      else
         call vberr('wrong call of diirw')
      end if
c
      return
      end
_ENDIF

c***********************************************************************
      subroutine diiseq(dm,lll,mmm,n,determ,dmt)
c
      implicit REAL (a-h,o-z), integer (i-n)
INCLUDE(common/scftvb)
INCLUDE(common/vbcri)
INCLUDE(../m4/common/iofile)
c
c...  set up diis equations and solve them
c...  results are returned in 1st n elements of dmt
c
      dimension dm(n,n),lll(n),mmm(n),dmt(*)
c
      if (iprins.gt.1) then
         write(iwr,11)
11       format(' DIIS - triangle')
         call tripr2(dmt,n)
      end if
      call square(dm,dmt,n,n)
c
c...  Now the matrix dm(n,n) contains <g(i)!g(j)> where i,j=1
c...  corresponds to the oldest value, 2 to the next one, etc.
c
c...  multiply diagonals by 1.01 (cf  JCP 84,5728 (1986)
c...  this indeed stabilises the process
c
      do 5 i=1,n
5     dm(i,i) = dm(i,i)*1.01d0
c
c...  Row i is subtracted from row i+1.  (compact equations)
c...  original equations are
c...  B11 B12 ... B1m -1       c1      0
c...  B21 B22 ... B2m -1       c2      0
c...  ... ... ... ... ..       ..      .
c...  Bm1 Bm2 ... Bmm -1       cm      0
c...   1   1       1   0      lamda    1
c...  (cf. P.Pulay CPL 73,393 (1980)) (note sign difference)
c...  substraction gets rid of langrange multplier lamda and
c...  reduces dimension by 1
c
      do 30 iii=n,2,-1
c...       substract
         do 10 i=1,n
10       dm(iii,i) = dm(iii,i) - dm(iii-1,i)
c...       normalise
         ss = dm(iii,iii)
         if (dabs(ss).lt.crilow) ss=crilow
         ss = 1.0d0/ss
         do 20 i=1,n
20       dm(iii,i) = dm(iii,i)*ss
30    continue
      do 40 i=1,n
40    dm(1,i) = 1.0d0
c
c...    now invert the diism_vb matrix. Its first column will give the coefficients
      call osinv_vb(dm,n,determ,crilow,lll,mmm)
c
c...    move result  to result => dmt
c
      do 50 i=1,n
50    dmt(i) = dm(i,1)
c
      return
      end

c***********************************************************************
      subroutine diism_vb(g,cr,ng,dm,nd,iu2,nskip,iu3)
c
      implicit REAL (a-h,o-z), integer (i-n)
c
c...  add to an old DIIS-matrix the extra row for the cuurent
c...  error-vector g (yielding a matrix of dimension nd)
c...  if nskip >0  we have to delete the oldest column
c...  iu2,iu3 : tape-numbers for old error-vectors/DIIS matrix
c...  nskip # records (old vectors/error-vectors) to skip
c
      dimension g(ng),cr(ng),dm(*)
INCLUDE(common/scra7vb)
INCLUDE(../m4/common/iofile)
c
      if ((iu2.ne.2.or.iu3.ne.3).and.nd.eq.1) 
     1   call caserr('diis file confusion')
c
c...  get old DIIS-matrix
c
      if (nd.gt.1) then
_IF(atmol)
         rewind iu3
         read(iu3)
         read(iu3)
_ENDIF
         if (nskip.le.0) then
_IF(atmol)
             call diirw(dm,nd*(nd-1)/2,iu3,'read')
_ELSE
             kvb3 = kscra7vb('kvb3',0,'r','n')
             call rdedx(dm,nd*(nd-1)/2,kvb3+2*nvb3,num8)
_ENDIF
         else
_IF(atmol)
             call diirw(dm,nd*(nd+1)/2,iu3,'read')
_ELSE
             kvb3 = kscra7vb('kvb3',0,'r','n')
             call rdedx(dm,nd*(nd+1)/2,kvb3+2*nvb3,num8)
_ENDIF
c...      delete oldest column
            kko = 0
            kkn = 0
            do 10 i=1,nd
              kko = kko + 1
              do 10 j=2,i
                kkn = kkn + 1
                kko = kko + 1
                dm(kkn) = dm(kko)
10          continue
         end if
      end if
c
c...  add new row
c
_IF(atmol)
      rewind iu2
      do 20 i=1,nskip
         read(iu2)
20    continue
_ENDIF
c
      kk = nd*(nd-1)/2
      do 30 i=1,nd-1
         kk = kk + 1
_IF(atmol)
         read(iu2) cr
_ELSE
         kkd = mod(i-1+nskip,maxdis)
         kvb2 = kscra7vb('kvb2',0,'r','n')
         call rdedx(cr,ng,kvb2+kkd*nvb12,num8)
_ENDIF
         dm(kk) = ddot(ng,cr,1,g,1)
30    continue
      dm(nd*(nd+1)/2) = ddot(ng,g,1,g,1)
c
c...  dump current diis matrix
c
_IF(atmol)
      rewind iu3
      read(iu3)
      read(iu3)
      call diirw(dm,nd*(nd+1)/2,iu3,'write')
_ELSE
      kvb3 = kscra7vb('kvb3',0,'r','n')
      call wrt3(dm,nd*(nd+1)/2,kvb3+2*nvb3,num8)
_ENDIF
c
c...  dump current (last) error-vector
c
_IF(atmol)
      write(iu2) g
_ELSE
      kkd = mod(nd-1+nskip,maxdis)
      kvb2 = kscra7vb('kvb2',0,'r','n')
      call wrt3(g,ng,kvb2+kkd*nvb12,num8)
_ENDIF
c
      return
      end

c***********************************************************************
      subroutine diisph(vo,vn,ndim,na,text)
c
      implicit REAL (a-h,o-z)
c
c...   check phases of vn vs. vo  (and overlaps)
c...   *note* all in orthogonalised basis
c...   do not normalise !!
c...   ** debug ** probably not needed
c
      dimension vo(ndim,na),vn(ndim,na)
INCLUDE(../m4/common/iofile)
      character*(*) text
c...  Phase should only be checked from second run of diis
      save icount
      data icount/0/
c
      icount=icount+1
      do 100 i=1,na
c
         ss = ddot(ndim,vn(1,i),1,vn(1,i),1)
         ss = 1.0d0/dsqrt(ss)
         call dscal(ndim,ss,vn(1,i),1)
c
         ss = ddot(ndim,vo(1,i),1,vn(1,i),1)
         if (ss.lt.0.0d0) then
            ll = len(text)
            if (icount.eq.1) then
              write(iwr,*) ' phase compare ',text(1:ll),
     *                 ' error at orbital ',i,' of ',ss,
     *                 'continuing because it is the first time'
            else
               write(iwr,*) ' phase compare ',text(1:ll),
     *                 ' error at orbital ',i,' of ',ss
               call vberr('phase error')
            end if
         else if (ss.lt.0.8d0) then
            ll = len(text)
            write(iwr,*) ' phase compare ',text(1:ll),
     *              ' ** warning overlap at mo ',i,' only ',ss
         end if
100   continue
c
      return
      end

c***********************************************************************
      subroutine diisrs(v,cr,n,c,nc,iun,nskip)
c
      implicit REAL (a-h,o-z), integer (i-n)
c
c...    calculate extrapolated c as linear combination of the c's
c...    v = sum v(i)*c(i)
c...   all vectors on unit iun
c...    nskip : # records to skip first
c
      dimension v(n),cr(n),c(nc)
INCLUDE(common/scra7vb)
_IFN(atmol)
INCLUDE(../m4/common/iofile)
      character*4 ckvb
c
      if (iun.eq.1) then
         ckvb = 'kvb1'
      else if (iun.eq.2) then
         ckvb = 'kvb2'
      else
         call caserr('diis file foulup')
      endif
_ELSE
c
      rewind iun
      do 10 i=1,nskip
10    read(iun)
_ENDIF
c
      call vclr(v,1,n)
c
      do 20 i=1,nc
_IF(atmol)
         read(iun) cr
_ELSE
         kvb = kscra7vb(ckvb,0,'r','n')
         kkd = mod(i-1+nskip,maxdis)*nvb12+kvb
         call  rdedx(cr,n,kkd,num8)
_ENDIF
         call daxpy(n,c(i),cr,1,v,1)
20    continue
c
      return
      end

c***********************************************************************
      subroutine sminh(s,t,n,val)
c
      implicit REAL (a-h,o-z), integer (i-n)
c
c...  determine the square s**-0.5 matrix t which transforms
c...  to a symmetrically orthogonal basis
c
      dimension s(*),t(n,n),val(n)
INCLUDE(common/vbcri)
c
      ntri = n*(n+1)/2
c
c     diagonalize the s-matrix
c
      call jaclow(s,t,n,ntri,crilow)
c
c     form s**(-.5)
c
c     make a row of the negative squareroot of the eigenvalues
c
      k = 0
      do 10 i=1,n
         k = k + i
         if (dabs(s(k)).lt.crilow) then
c....       eigenvector effectively zero
            val(i) = 0.0d0
            call vberr('dependent basis in sminh')
         else
            val(i) = 1/dsqrt(s(k))
         end if
10    continue
c
c     form t.val.t-cross  (triangle)
c
      call vclr(s,1,ntri)
      do 30 k=1,n
         m=0
         do 20 i=1,n
            call daxpy(i,t(i,k)*val(k),t(1,k),1,s(m+1),1)
            m=m+i
20       continue
30    continue
c
c      square
c
      call square(t,s,n,n)
c
      return
      end

c***********************************************************************
      subroutine osinv_vb (a,n,d,tol,l,m)
      implicit REAL (a-h,o-z)
c
c     parameters:  a - input matrix , destroyed in computation and repla
c                      by resultant inverse (must be a general matrix)
c                  n - order of matrix a
c                  d - resultant determinant
c            l and m - work vectors of lenght n
c                tol - if pivot element is less than this parameter the
c                      matrix is taken for singular (usually = 1.0e-8)
c     a determinant of zero indicates that the matrix is singular
c
      dimension a(*), m(*), l(*)
      d=1.0d0
      nk=-n
      do 180 k=1,n
         nk=nk+n
         l(k)=k
         m(k)=k
         kk=nk+k
         biga=a(kk)
         do 20 j=k,n
            iz=n*(j-1)
         do 20 i=k,n
            ij=iz+i
c
c     10 follows
c
            if (dabs(biga)-dabs(a(ij))) 10,20,20
   10       biga=a(ij)
            l(k)=i
            m(k)=j
   20    continue
         j=l(k)
         if (j-k) 50,50,30
   30    ki=k-n
         do 40 i=1,n
            ki=ki+n
            holo=-a(ki)
            ji=ki-k+j
            a(ki)=a(ji)
   40    a(ji)=holo
   50    i=m(k)
         if (i-k) 80,80,60
   60    jp=n*(i-1)
         do 70 j=1,n
            jk=nk+j
            ji=jp+j
            holo=-a(jk)
            a(jk)=a(ji)
   70    a(ji)=holo
   80    if (dabs(biga)-tol) 90,100,100
   90    d=0.0d0
         return
  100    do 120 i=1,n
            if (i-k) 110,120,110
  110       ik=nk+i
            a(ik)=a(ik)/(-biga)
  120    continue
         do 150 i=1,n
            ik=nk+i
            ij=i-n
         do 150 j=1,n
            ij=ij+n
            if (i-k) 130,150,130
  130       if (j-k) 140,150,140
  140       kj=ij-i+k
            a(ij)=a(ik)*a(kj)+a(ij)
  150    continue
         kj=k-n
         do 170 j=1,n
            kj=kj+n
            if (j-k) 160,170,160
  160       a(kj)=a(kj)/biga
  170    continue
         d=d*biga
         a(kk)=1.0d0/biga
  180 continue
      k=n
  190 k=k-1
      if (k) 260,260,200
  200 i=l(k)
      if (i-k) 230,230,210
  210 jq=n*(k-1)
      jr=n*(i-1)
      do 220 j=1,n
         jk=jq+j
         holo=a(jk)
         ji=jr+j
         a(jk)=-a(ji)
  220 a(ji)=holo
  230 j=m(k)
      if (j-k) 190,190,240
  240 ki=k-n
      do 250 i=1,n
         ki=ki+n
         holo=a(ki)
         ji=ki+j-k
         a(ki)=-a(ji)
  250 a(ji)=holo
      go to 190
  260 return
c
      end

c***********************************************************************
      subroutine fmos(smo,sao,c,v,nmdim,ndim,epsss)
c
      implicit REAL (a-h,o-z), integer (i-n)
INCLUDE(../m4/common/iofile)
c
c...  compute the s-matrix on mo-basis for nmdim mo's (c)
c...  of dimension ndim with ao-overlap-matrix sao
c
      dimension smo(1),sao(1),c(1),v(1)
c
      ks = 1
      ki = 1
c
      do 20 i=1,nmdim
         call cntrc(sao,c(ki),v,ndim,epsss)
         kj = 1
         do 10 j=1,i
            smo(ks) = ddot(ndim,c(kj),1,v,1)
            ks = ks + 1
10       kj = kj + ndim
20    ki = ki + ndim
c
      return
      end

c***********************************************************************
      integer function iflip(ifrom,ito)
      implicit REAL  (a-h,o-z)  , integer  (i-n)
c.....
c.....this returns the sign associated with the brillouin excitation
c.....ifrom => ito, as defined in the input of the mix-option
c.....using this facility e.g. one can make two sp2-like hybrids in
c.....h2o (pointing at the h-s) equivalent (they use one p-orbital
c.....equally, but with the opposite sign)
c.....
INCLUDE(../m4/common/iofile)
INCLUDE(../m4/common/sizes)
INCLUDE(common/turtleparam)
INCLUDE(common/brill)
INCLUDE(common/vbequiv)
      common /fast/ nospin,nosymm,nosymc
      logical nospin,nosymm,nosymc
c
      save iagain
      data iagain/0/
      if (equiv.and.iagain.eq.0) write(iwr,'(a)') 
     1  ' **** iflip not understood for hand equivalence ****'
      iagain = 1
c.....locate ifrom in iexw
      if (ito.le.nscf.and.1.eq.0) then
c..
c.....   this is an internal excitation, no sign-flips (taken care of by
c.....   "phase")
c..
         iflip = 1
         return
      end if
      ipos1 = 0
      ipos2 = 0
      do 10 i=1,nscf
         if (iex(1,i).eq.ifrom) ipos1 = i
10    continue
      if (ipos1.eq.0) then
         write(iwr,123)ifrom,ito
123      format(1x,i3,' cannot be found on mo list. ito : ',i3)
         call vberr('? iflip/1 ?')
      end if
c.....locate ito
      do 30 i=2,nex(ipos1)+1
         if (iex(i,ipos1).eq.ito) ipos2 = i
30    continue
      if (ipos2.eq.0) then
         write(iwr,234) ito,ifrom
234      format(1x,i3,' cannot be found in ',i3,1h','s excitations,')
         call vberr('? iflip/2 ?')
      end if
      iflip = ieqsig(ipos2,ipos1)
c...    iflip is unclear for equic
      if (equiv) iflip = 1
      return
      end

c***********************************************************************
      integer function ioccup(iorb,iconf,nelec)
c.....
c.....return the occupation number of iorb in iconf
c.....
      implicit REAL (a-h,o-z)
      dimension iconf(nelec)
      icount = 0
      do 10 i=1,nelec
         if (iconf(i).eq.iorb) icount = icount + 1
10    continue
      ioccup = icount
      return
      end

c***********************************************************************
      subroutine jaccov(tens,n,s,vec,e,iorder,crit,scr)
      implicit REAL (a-h,o-z), integer(i-n)
c.....
c.....find eigenvectors of a covariant first order
c.....hermitian (triangular storage) tensor tens, with metric s
c.....tens : tensor
c.....n : dimension
c.....s : metric  ( will be spoiled !!!! )
c.....vec : eigenvectors
c.....e : eigenvalues
c.....iorder : order-parameter for jacobi (see that routine)
c.....crit   : criterion        ,,  ,,           ,,
c.....scr : scratch
c.....
      dimension tens(*),s(*),vec(n,n),e(*),scr(*)
      call vclr(vec,1,n*n)
      do i=1,n
         vec(i,i)=1.0d0
      enddo
      kscr = n*(n+1)/2 + 1
      call trtrtr(s,tens,scr,n,scr(kscr))
      kto = kscr
      ktrd = kto + n*(n+1)/2
      kscr = ktrd + n*n
      call schmidt(s,vec,n,itroub)
      if (itroub.ne.0) then
         call tripri(s,n)
         call vberr('rotten metric in jaccov')
      end if
      call mult11(scr,scr(kto),scr(kscr),n,n,vec,scr(ktrd))
      kiky = ktrd
      call filiky(scr(kiky),n+1)
      call jacobt(scr(kto),scr(kiky),n,vec,n,e,1,iorder,crit,
     &                                                 scr(kscr))
      return
      end

c***********************************************************************
      subroutine jaclow(a,b,n,nnhalf,eps)
c
      implicit REAL  (a-h,o-z) , integer   (i-n)
c
c
c..   jacobi diagonaliser for lowdin orthogonalisation
c
      dimension a(nnhalf), b(n,n)
      REAL ii,jj,ij
      integer contr
      nmin1 = n - 1
      do 10 i=1,n
      do 10 j=1,n
10    b(i,j) = 0.d0
      do 20 i=1,n
20    b(i,i) = 1.d0
      if (n.le.1) return
100   u = 0
      l = 2
      do 110 i = 2,n
      do 120 j = 2,i
      u = u + a(l)**2
      l = l+1
120   continue
      l = l+1
110   continue
      u = dsqrt(u/nnhalf)
      if (u.lt.eps) goto 300
      i = 1
      j = 1
      kii = 1
      kjj = 1
      kij = 1
      contr = 0
      goto 200
202   if (j.eq.nmin1) goto 400
      i = j+1
      j = i
      kii = kjj + j
      kjj = kii
      kij = kii
200   if (i.eq.n) go to 202
      kij = kij + i
      i = i + 1
      kii = kii + i
      ij = a(kij)
      if (dabs(ij).lt.u) goto 200
      contr = contr + 1
      ii = a(kii)
      jj = a(kjj)
      vm = (jj-ii)/2.d0
      vn = dsqrt(ij*ij+vm*vm)
      cos = (vn+dabs(vm))/(2.d0*vn)
      co = dsqrt(cos)
      si = ij/(dsign(2.d0,vm)*vn*co)
      sin = si*si
      imin = i-1
      jmin = j-1
      iplus = i + 1
      jplus = j+1
      ik = kii - imin
      kj = kjj - jmin
      if (j.eq.1) goto 211
      do 210 kv = 1,jmin
      v = a(ik)
      w = a(kj)
      a(ik) = co*v - si*w
      a(kj) = si*v + co*w
      ik = ik+1
      kj = kj+1
210   continue
211   akjj = ii * sin + jj * cos + 2 * si * co * ij
      a(kjj) = akjj
      a(kij) = 0
      ik = ik + 1
      kj = kj + j
      if(j.eq.imin) goto 213
      do 212 kv = jplus , imin
      v = a(ik)
      w = a(kj)
      a(ik) = co*v - si*w
      a(kj) = si*v + co*w
      ik = ik+1
      kj = kj+kv
212   continue
213   a(kii) = ii + jj - akjj
      ik = ik + i
      kj = kj + i
      if(i.eq.n) goto 215
      do 214 kv = iplus , n
      v = a(ik)
      w = a(kj)
      a(ik) = co*v - si*w
      a(kj) = si*v + co*w
      ik = ik+kv
      kj = kj+kv
214   continue
215   do 216 kv = 1,n
      v = b(kv,i)
      w = b(kv,j)
      b(kv,i) = co*v - si*w
      b(kv,j) = si*v + co*w
216   continue
      goto 200
400   u = u / 6
      if(contr.lt.4) u = u/2
      if (u.lt.eps) go to 100
      contr = 0
      i = 1
      j = 1
      kii = 1
      kjj = 1
      kij = 1
      goto 200
300   return
      end

c***********************************************************************
      subroutine jacobt(h,iky,nbasis,v,nrow,e,init,iorder,small,y)
c
      implicit REAL  (a-h,o-z) , integer   (i-n)
c
c.....
c.....standard jacobi diagonaliser. h is diagonalised, v will contain
c.....the eigenvectors, e the eigenvalues. iky should contain i*(i-1)/2
c.....at i. y is a scratch array of length 2*nbasis.
c.....init   : = 1 ? => v is assumed to be initialised elsewhere
c.....  ,,   : other values will cause it to be a unit matrix at first
c.....iorder : = 0 ? => no ordering in the vectors, no e.val. at all
c.....  ,,   : = 1 ? => no ordering in the vectors, e.val's in e.
c.....  ,,   : = 2 ? => increasing eigenvalues/vectors
c.....  ,,   : = 3 ? => decreasing eigenvalues/vectors
c.....  ,,   : = 4 ? => lock mode (i.e. try to change the initial
c.....                             vectors the least possible)
c.....
      dimension h(*),iky(*),v(nrow,nbasis),e(*),y(*)
_IF(parallel)
c..
c.. Was meant to have results the same for all nodes. This is now
c.. done by ordering the vectors not only to eigenvalue, but
c.. also by vector elements. If problems with different vectors
c.. on different nodes occur this can be reenabled. If it is
c.. not used less barriers. (Fokke 2001)
c      logical oroot
_ENDIF
      if (init.ne.1) then
         call vclr(v,1,nrow*nbasis)
         do 10 i=1,nbasis
            v(i,i) = 1.0d0
10       continue
      end if
      if (nbasis.eq.1) then
         e(1)   = h(1)
         goto 999
      end if
_IF(parallel)
c     if (oroot()) then
_ENDIF
11    rlarge = 0.0d0
      do 30 i=2,nbasis
         ii = iky(i)
         do 20 j=1,i-1
           if (dabs(h(ii+j)).gt.rlarge) rlarge = dabs(h(ii+j))
20       continue
30    continue
      if (rlarge.gt.small) then
         reason = rlarge * .1d0
         do 60 i=1,nbasis - 1
            call fmove( h(iky(i)+1),y,i)
            call gather(nbasis-i,y(i+1),h(i+1),iky(i+1))
            hii = y(i)
            do 50 j=i+1,nbasis
               if (dabs(y(j)).ge.reason) then
                  vii  = hii * 0.5d0
                  vij  = y(j)
                  vjj  = h(iky(j)+j) * 0.5d0
                  diff = vii - vjj
                  root = dsqrt( diff*diff + vij*vij )
                  if (diff.lt.0.0d0) root = -root
                  cosr = (diff + root)/vij
                  sinr = dsqrt( 1.0 d0/ (1.0d0+cosr*cosr) )
                  cosr = cosr * sinr
                  diff = vii + vjj
                  call drot(j-1,y,1,h(iky(j)+1),1,cosr,sinr)
                  hii  = diff + root
                  h( iky(j)+j ) = diff - root
                  y(j) = 0.0d0
                  do 40 k=j+1,nbasis
                     vij   = y(k)
                     kj    = iky(k) + j
                     y(k)  = vij   * cosr + h(kj) * sinr
                     h(kj) = h(kj) * cosr - vij   * sinr
40                continue
                  call drot(nbasis,v(1,i),1,v(1,j),1,cosr,sinr)
               end if
50          continue
            call fmove(y,h(iky(i)+1),i-1)
            call scatter(nbasis-i,h(i+1),iky(i+1),y(i+1))
            h(iky(i)+i) = hii
60       continue
         goto 11
      end if
      if (iorder.eq.0) then
         goto 999
      else if (iorder.eq.1) then
         call gather(nbasis,e,h,iky(2))
         goto 999
      else if (iorder.eq.2) then
         call gather(nbasis,e,h,iky(2))
         call ordert(e,v,nbasis,nrow,y,y(nbasis+1), 1)
         call scatter(nbasis,h,iky(2),e)
      else if (iorder.eq.3) then
         call gather(nbasis,e,h,iky(2))
         call ordert(e,v,nbasis,nrow,y,y(nbasis+1),-1)
         call scatter(nbasis,h,iky(2),e)
      end if
999   continue
_IF(parallel)
c     endif
c     call pg_brdcst(100,h,nbasis*(nbasis+1)/2,0)
c     if (iorder.ne.0) call pg_brdcst(200,e,nbasis,0)
c     call pg_brdcst(300,v,nbasis*nrow,0)
_ENDIF
      return
      end

c***********************************************************************
      subroutine ordert(e,v,n,nrow,y,imap,ipar)
c
      implicit REAL  (a-h,o-z) , integer   (i-n)
c
      dimension e(*),v(nrow,n),y(*),imap(*)
      do 10 i=1,n
         imap(i) = i
10    continue
      call bubblevec(e,n,v,nrow,ipar,imap)
      do 20 iold=1,n-1
         inew = imap(iold)
         imap(iold) = 0
         if (iold.eq.inew) goto 20
         call fmove(v(1,iold),y,n)
         call fmove(v(1,inew),v(1,iold),n)
         call fmove(y,v(1,inew),n)
         imap(locati(iold,imap,n)) = inew
20    continue
      return
      end

c***********************************************************************
      subroutine jacobs(h,n,s,vec,e,iorder,crit,scr)
      implicit REAL (a-h,o-z), integer(i-n)
c.....
c.....solve generalised eigenvalue problem,h is spoiled (by mult11)
c.....h : triangular matrix
c.....n : dimension
c.....s : metric (triangle), is spoiled by schmidt
c.....vec : eigenvectors
c.....e : eigenvalues
c.....iorder,crit : parameters for jacobi (see that routine)
c.....scr : scratch
c.....
      dimension h(*),s(*),vec(n,n),e(*),scr(*)
      kto = 1
      ktrd = kto + n*(n+1)/2
      kscr = ktrd + n*n
      call schmidt(s,vec,n,itroub)
      if (itroub.ne.0) call vberr('rotten metric in jacobs')
      call mult11(h,scr(kto),scr(kscr),n,n,vec,scr(ktrd))
      kiky = ktrd
      call filiky(scr(kiky),n+1)
      call jacobt(scr(kto),scr(kiky),n,vec,n,e,1,iorder,crit,
     &                                                 scr(kscr))
      return
      end
 
c***********************************************************************
      subroutine normt(v,nbasis,ncol,s,w)
c
      implicit REAL  (a-h,o-z) , integer   (i-n)
INCLUDE(../m4/common/sizes)
INCLUDE(common/vbcri)
c
c...   normalise ncol vectors of dimension nbasis with s-matrix s
c...   allow zero vectors but have to be last and reduce dimension 
c
      dimension v(nbasis,ncol),s(*),w(nbasis)
c
      nccc = ncol
c
      do 100 i=1,ncol
         call cntrc(s,v(1,i),w,nbasis,crilow)
         ss = ddot(nbasis,v(1,i),1,w,1)
         if (ss.lt.crilow) then
            nccc = nccc - 1
         else
            if (nccc.ne.ncol) 
     1      call caserr('harmon confusion; 0-vectors not last')
            tnorm = 1.0d0/dsqrt(ss)
            do j=1,nbasis
               v(j,i) = v(j,i)*tnorm
            end do
         end if
100   continue
c
      ncol = nccc
c
      return
      end

c***********************************************************************
      subroutine oldnew(bold,bnew,ndet,idet,nelec,nalfa,
     &                  vold,vnew,sao,smo,scr,nbasis,nact,swave)
      implicit REAL  (a-h,o-z) , integer   (i-n)
c.....
c.....determine the overlap of the old wave-function with the new one
c.....
      dimension bold(ndet),bnew(ndet),idet(nelec,ndet),
     &          vold(nbasis,nact),vnew(nbasis,nact),
     &          sao(nbasis,nbasis),smo(nact,nact),scr(*)
c.....
c.....symmetrise s
c.....
      kscr2 = nbasis*nact+1
      call trisqu(sao,scr(kscr2),nbasis)
c.....make non-hermitian overlap-matrix between old and new orbitals
c.....
      call vclr(scr,1,nbasis*nact)
      call mxmb(vold,nbasis,1,scr(kscr2),1,nbasis,scr,1,nact,
     &                                  nact,nbasis,nbasis)
      call vclr(smo,1,nact*nact)
      call mxmb(scr,1,nact,vnew,1,nbasis,smo,1,nact,
     &                                  nact,nbasis,nact)
c.....
c.....smo now contains the overlap-matrix between the old and the
c.....new orbitals
c.....
      ibold = 0
      swave = 0.0d0
      do 70 i=1,ndet
         ibold = ibold + 1
         ibnew = 0
         do 60 j=1,ndet
            ibnew = ibnew + 1
            it = 1
            do 30 k=1,nalfa
               do 20 l=1,nalfa
                  kk = idet(k,i)
                  ll = idet(l,j)
                  scr(it) = smo(kk,ll)
                  it = it + 1
20             continue
30          continue
            call determ(scr,nalfa,deta,irank)
            it = 1
            do 50 k=nalfa+1,nelec
               do 40 l=nalfa+1,nelec
                  kk = idet(k,i)
                  ll = idet(l,j)
                  scr(it) = smo(kk,ll)
                  it = it + 1
40             continue
50          continue
            call determ(scr,nelec-nalfa,detb,irank)
            swave = swave + deta*detb*bold(ibold)*bnew(ibnew)
60      continue
70    continue
      swave = dabs(swave)
      return
      end

c***********************************************************************
      subroutine orbopt(v,nb,nact,vnew,bi,q,lword)
      implicit REAL (a-h,o-z), integer (i-n)
c
c...  use info in ci-vector to update orbitals
c...  if converged return scfconv = .true.
c...  v are all the orbitals (including core)
c...  bi ci-vector
c...  q,lq (same) scratch space
c
      dimension v(nb,*),bi(*),q(lword),vnew(nb,nact)
      character*10 charwall
c
INCLUDE(common/vbdiist)
INCLUDE(common/tractlt)
INCLUDE(common/hsinfo)
c
c...  nstruc (in hsinfo) is now accidently dimension of bi(*)
c
INCLUDE(common/turtleparam)
INCLUDE(../m4/common/sizes)
INCLUDE(common/vbpert)
INCLUDE(../m4/common/iofile)
INCLUDE(common/brill)
INCLUDE(common/scftvb)
INCLUDE(common/vbcri)
INCLUDE(common/splice)
INCLUDE(common/infato)
INCLUDE(common/vbequiv)
      dlog10(aa) = dlog(aa) / dlog(10.d0)
c
c...  check convergence
c
c...  make C0 > 0 to avoid phase flips (diis)
c
      scfconv = .false.
      div = 0.0d0
      do 10 i=2,nstruc
         if (bi(1).lt.0.0d0) bi(i) = -bi(i)
10    div = dmax1(dabs(bi(i)),div)
      if (bi(1).lt.0.0d0) bi(1) = -bi(1)
      if (optcri.eq.2) then
c.....   use overlap of current wavefunction with previous one as crit.
         if (swave.ge.criscf) scfconv = .true.
      else if (optcri.eq.1) then
c.....   use largest change in mo coefficient as criterion
         if (div.le.criscf) scfconv = .true.
      else 
c.....   use real brillouin theorem
         if (brmax.le.criscf) scfconv = .true.
      end if
      if (nitscf.ge.maxscf) then
         scfconv = .true.
         write(iwr,604) nitscf
604   format(/,' iterations ended due to excessive iterations, viz.',i4)
      end if
      del = dabs(1.d0-swave)
c.....1.0e-14 is machine precision .. we may take critor
      if (del.gt.critor) then
         del = dlog10(del) - .5d0
         idel=idint(del)
      else
         idel = -15
      end if
      if (nitdiis.eq.0) then
        if (ovbper) then
          write(iwr,611) nitscf,cpulft(1),charwall(),evb,ebi+shiscf,
     1                   div,idel,brmax
611       format(/,' it.',i5,' at',f9.2,a10,' evb',f13.7,
     &             ' ebi-p',f13.7,' div ',1pe7.1,'/del',i3,'/brm '
     &               ,1pe7.1)
        else
          write(iwr,600) nitscf,cpulft(1),charwall(),evb,ebi+shiscf,
     1                   div,idel,brmax
600       format(/,' it.',i5,' at',f9.2,a10,' evb',f13.7,' ebi',f13.7,
     &             ' div ',1pe7.1,'/del',i3,'/brm ',1pe7.1)
        end if
      else
        if (ovbper) then
          write(iwr,622) nitscf,cpulft(1),charwall(),evb,
     &                   ebi+shiscf,div,idel,
     &            brmax,nitdiis,int(perccur),maxdiis,int(percmax)
622       format(/,' it.',i5,' at',f9.2,a10,' evb',f13.7,
     &             ' ebi-p',f13.7,' div ',1pe7.1,'/del',i3,'/brm ',
     &        1pe7.1,'  diis',i3,' =',i3,' %',i2,' =',i3,'%')
        else
          write(iwr,633) nitscf,cpulft(1),charwall(),evb,ebi+shiscf,
     +           div,idel,brmax,nitdiis,int(perccur),maxdiis,
     +           int(percmax)
633       format(/,' it.',i5,' at',f9.2,a10,' evb',f13.7,' ebi',f13.7,
     &             ' div ',1pe7.1,'/del',i3,'/brm ',1pe7.1,
     &      '  diis',i3,' =',i3,' %',i2,' =',i3,'%')
        end if
      end if
c
c...  build transformation matrix
c
      call vclr(q,1,nsa*nact+nsa*nact)
      do 6502 i=1,nact
6502  q((i-1)*nsa+i) = 1.0d0
c
      iq = 1
      jq = 1
      iibi = 1
      maxoc = 0
      do 44 i=1,nscf
         ifrom = iex(1,i)
c
c...  if ifrom is 0 it means that this orbital is frozen so skip it
c
         if (ifrom.ne.0) then
            maxoc = max(maxoc,ifrom)
            q( (ifrom-1)*nsa + ifrom ) = bi(1)
            if (fixb0.ge.0.0d0) q( (ifrom-1)*nsa + ifrom ) = fixb0
            do 33 j=2,nex(i)+1
               iibi = iibi + 1
               ito = iex(j,i)
               q( (ifrom-1)*nsa + ito ) = bi(iibi) * iflip(ifrom,ito) 
     1                                             * dampvb
33          continue
            jq = jq + 1
            if (jq.le.iequi(iq)) then
               iibi = iibi - nex(i)
            else
               iq = iq + 1
               jq = 1
            end if
         end if
44    continue
c
c...  update orbitals
c
      if (equiv) call primequiv(v,nbasis,1,'occ','beg orbopt')
      if (iprins.gt.20) then
         write(iwr,11)
11       format(/,' orbitals before update, mos =>  ')
         call prvc(v(1,ncore+1),nact,nbasis,v,'o','l')
         write(iwr,22)
22       format(/,' transformation matrix  ')
         call prvc(q,maxoc,nsa,q,'o','a')
      end if
c
c...  dump old vectors to tape33 for oldnew (determines overlap between
c...  old and new wavefunctions) and natorb (1-electron density matrix)
c
_IF(atmol)
      write(33) nb,ncore+nact,((v(ii,jj),ii=1,nb),jj=1,ncore+nact)
_ELSE
      kvb7_vo = kscra7vb('kvb7_vo',nb*(ncore+nact),'r','w')
      kvb7_vn = kscra7vb('kvb7_vn',nb*(ncore+nact),'r','w')
      kvb7_vnew = kscra7vb('kvb7_vnew',nb*nact,'r','w')
      call wrt3(v,nb*(ncore+nact),kvb7_vo,num8)
_ENDIF
c
c...  if maxscf=0 we have done the work but do not change the orbitals
c
      if (maxscf.eq.0) then
         write(iwr,23)
23       format(' ****** Orbitals are unchanged ******')
         call fmove(v(1,ncore+1),vnew,nbasis*nact)
      else
       call mxmd(v(1,ncore+1),1,nbasis,q,1,nsa,vnew,1,nbasis,
     &           nbasis,nsa,nact)
      end if
c
c...  pert (nowadays opti)
c
      if (idel.le.idelp) ovbper =.true.
c
c...  diis (normalises q itself)
c
      if (idel.lt.iwdiis.and.nitscf.ge.ntdiis.and..not.scfconv) then
         ovbdiis = .true.
         kscr = 2*nsa*nact+1
         call vbdiis(q,v(1,ncore+1),vnew,nbasis,nact,
     &               q(kscr),lword-kscr,madiis,midiis)
c
c...     the diis transformations may ruin hybrids, so re-clear
c
         if (clean.and.hybry) then
            call clvecvb(v(1,ncore+1),nact,nbasis,'active','small',
     1                   'after diis')
         end if
      else
         call fmove(vnew,v(1,ncore+1),nbasis*nact)
      end if
c
c...  normalise orbitals and vnew as well / read s into q
c
      call get1e(q,dummy,'s',q)
      call normt(v(1,ncore+1),nbasis,nact,q,q(nbasis*(nbasis+1)/2+1))
      call normt(vnew,nbasis,nact,q,q(nbasis*(nbasis+1)/2+1))
c
c...  dump new vectors to tape33 for oldnew (determines overlap between
c...  old and new wavefunctions) and natorb (1-electron density matrix)
c...  info : orbitals updated but singles not orthogonal to doubles yet
c...  (caused bug in overlap old/new and natural orbitals) now in vbscf
c...  also the previous vnew vectors in case diis fails
c
      kvb7_vnew = kscra7vb('kvb7_vnew',nb*nact,'r','w')
      call wrt3(vnew,nb*nact,kvb7_vnew,num8)
c...     make sure equivalence restrictions are obeyed 
      if (equiv) call primequiv(v,nbasis,1,'occ','end orbopt')
c
_IF(atmol)
      call vberr(' vnew not properly  dumped')
      write(33) nb,ncore+nact,((v(ii,jj),ii=1,nb),jj=1,ncore+nact)
_ELSE
      kvb7_vn = kscra7vb('kvb7_vn',nb*(ncore+nact),'r','w')
      call wrt3(v,nb*(ncore+nact),kvb7_vn,num8)
_ENDIF
c
      if (iprins.ge.10) then
         write(iwr,602) (bi(i),i=1,nstruc)
602      format(/,' brillouin-state vector:',/,(10e12.4))
         write(iwr,601)
601      format(//,' --- new vb orbitals ---')
         call prvc(v(1,ncore+1),nact,nbasis,v,'o','l')
      end if
c
      return
      end

c***********************************************************************
      subroutine redorb(s,v,vcopy,vcore,iocvi,iset,iact,idoc,isoc,ivir,
     &                  hmoao,iex2,nbasis,qq)
      implicit REAL (a-h,o-z), integer (i-n)
c
c     v - active + virtual orbitals (nactiv+nbasis)
c     vcore - core orbitals (really all, as the core are just in front)
c 
INCLUDE(../m4/common/sizes)
INCLUDE(../m4/common/iofile)
INCLUDE(common/turtleparam)
c.....
      dimension s(nbasis*(nbasis+1)/2),v(nbasis,*),vcopy(nbasis,*),
     &          vcore(nbasis,*),iocvi(*),
     &          iset(nbasis),iact(nbasis),hmoao(*),
     &          iex2(nbasis+1,nscf),idoc(*),isoc(*),ivir(*),
     &          qq(*)
INCLUDE(common/scftvb)
INCLUDE(common/twice)
INCLUDE(common/splice)
INCLUDE(common/brill)
      logical there1,there2,crash
INCLUDE(common/infato)
INCLUDE(../m4/common/harmon)
INCLUDE(../m4/common/gmempara) 
      common /forbex/ nforbi,iforbi(maxact*maxact,2),
     &                nionj,ionj(maxact,2),
     &                nsetor,isetor(maxact+1,maxact),
     &                nionres,ionrest(maxact)
      common /fast/ nospin,nosymm,nosymc
      logical nospin,nosymm,nosymc
      dimension idref(mxorbvb)
      crash = .false.
c.....establish excitation patterns
      nao=nbasis
      nmo=nscf
      nactiv = nsingly + ndoubly
c
      call toetoe(v,hmoao,nao,nmo)
c.....
c.....   temporarily strip the metric (remove inter-atomic overlaps)
c.....   to prevent another center's virtual to resemble an atom's
c.....   occupied orbital
c.....
      if (hybry) call clmetr(s,nbasis,'any')
c.....
c.....iexw contains the original excitation-patterns, and is not changed
c.....iexw and nexw are set in toetoe
c.....iex will contain the adapted (renumbered,reduced) excitation
c.....patterns, that can be different in each iteration
c.....
      do 10 i=1,nscf
         do k=1,nexw(i)+1
            iex(k,i) = iexw(k,i)
         end do
         nex(i) = nexw(i)
10    continue
c
c.....iocvi contains per occupied the most resembling virtual
c.....order : frozen core, doubly occ. ,singly/variable occ.
c..... for the occupied orbitals, the order is given by idoubly,isingly
c..... as the core orbitals are not included in v, this gets weird
c.....Initially the virtuals span the complete functionspace.
c.....In the course of this routine the dimension of the virtual
c.....space will be diminished as much as possible.
c.....
      call izero(ncore+nactiv,iocvi,1)
c....
c.... find the virtuals that resemble the active orbitals
c.... *note* vcore is really core + activ
c....
      call fmove(v(1,nactiv+1),vcopy,nbasis*nbasis)
c.... vcopy are virtuals only (includes occupieds)
      kww  = igmem_alloc_inf(nbasis,'vbscf.m','redorb',
     @                       'kww',IGMEM_DEBUG)
      do i=1,nactiv+ncore
         iorb = iococ(i,'ord')
         maxs = mvocvi(vcore(1,iorb),vcopy,nbasis,nbasis,s,qq(kww))
         iocvi(i) = maxs + nactiv
c....       exclude the chosen virtual from the list of candidates
         call vclr(vcopy(1,maxs),1,nbasis)
      end do
      call gmem_free_inf(kww,'vbscf.m','redorb','kww')
      call fmove(v(1,nactiv+1),vcopy,nbasis*nbasis)
      do i=1,nmo
         ii = iococ(i+ncore,'ord') + nmo 
         if (locati(ii,iocvi(ncore+1),nmo).eq.0) call caserr(' misjoop')
      end do
c
c
c.....analyse the function-space
c
c      call reexc(hmoao,nao,nmo,v,iocvi,crash)
c
      call get1e(s,dummy,'s',s)
c      if (hybry) call get1e(s,dummy,'s',s)
      if (hybry) call clmetr(s,nbasis,'any')
c.....
      if (iprins.ge.100) then
         write(iwr,1) nactiv
1        format(/,' correspondence core-doubly-singly with virtuals ',
     1            '(shifted by',i4,')')
         if (ncore.gt.0) write(iwr,2)'core',(iococ(ij,'ord'),ij=1,ncore)
         if (ncore.gt.0) write(iwr,2)'virt',(iocvi(ij),ij=1,ncore)
         write(iwr,2) 'doub',(iococ(ij,'ord'),ij=ncore+1,ncore+ndoubly)
         write(iwr,2) 'virt',(iocvi(ij+ncore),ij=1,ndoubly)
	 write(iwr,2) 'sing',(iococ(ij+ncore+ndoubly,'ord'),ij=1,nsingly)
         write(iwr,2) 'virt',(iocvi(ij+ncore+ndoubly),ij=1,nsingly)
2        format(1x,a,(t6,40i3))
      end if
c.....if no orthogonalisations are requested skip that part
      if (.not.ortscf) goto 240
c.....
c.....analyse the subset structure of the excitation patterns
c.....perform all orthogonalisations that are allowed for, and then
c.....reduce the dimension of the orbitals space as much as possible
c.....
      if (iprins.gt.1000) write(iwr,11)
11    format(/,' subset analysis of the total orbital space :')
c.....
c.....copy iexw, remove doubly occupied virtuals, and change excitations
c.....to singly occupied virtuals into singly occupied occupieds
c.....(if .not.super)
c.....
      do 70 j=1,nscf
         iex2(1,j) = iexw(1,j)
         do 60 i=2,nexw(j)+1
            iex2(i,j) = iexw(i,j)
            do  k=1,ncore+ndoubly
               if (iex2(i,j).eq.iocvi(k)) iex2(i,j) = 0
            end do
            if (.not.super) then
c.....         exploit singly occupieds  (do not do this in case
c.....         orbitals are optimised per csf => super=.true.)
               do 50 k=1,nsingly
                  if (iex2(i,j).eq.iocvi(k+ncore+ndoubly)) 
     1                iex2(i,j) = isingly(k)
50             continue
            end if
60       continue
70    continue
      number = 0
      do 230 i=1,nscf
         do 220 j=2,nexw(i)+1
            iorb1 = iex2(j,i)
            if (iorb1.eq.0) goto 220
            number = number + 1
c.....      find subset belonging to iorb1
            call subex(iorb1,idoc,ndoc,isoc,nsoc,ivir,nvir,iex2,nbasis)
c.....
c.....      gather the doubly occupieds, singly occupieds and virtuals
c.....
            call vecgat(v,vcopy(1,1          ),nbasis,idoc,ndoc)
            call vecgat(v,vcopy(1,ndoc+1     ),nbasis,isoc,nsoc)
            call vecgat(v,vcopy(1,ndoc+nsoc+1),nbasis,ivir,nvir)
c.....      do orthogonalisations/projections
            call proorb(vcopy,s,nbasis,ndoc,nsoc)
            if (iprins.gt.1000000) then
c.....         check the result of proorb
               write(iwr,22) number
22             format(' ---------- set # ',i2,' ----------')
               write(iwr,33) (idoc(ij),ij=1,ndoc)
               write(iwr,44) (isoc(ij),ij=1,nsoc)
               write(iwr,55) (ivir(ij),ij=1,nvir)
33             format(/,' doubly occupied  :',(t20,30i3))
44             format(/,' singly occupied  :',(t20,30i3))
55             format(/,' virtual orbitals :',(t20,30i3))
               write(iwr,'(a)') ' checking orthogonality doubly' 
               call ortche(vcopy,ndoc,
     &                     vcopy,ndoc,nbasis,s,0,qq)
               write(iwr,'(a)') ' checking orthogonality singly' 
               call ortche(vcopy(1,ndoc+1),nsoc,
     &                     vcopy,nsoc,nbasis,s,0,qq)
               write(iwr,'(a)') ' checking orthogonality virtuals' 
               call ortche(vcopy(1,ndoc+nsoc+1),nvir,
     &                     vcopy,nvir,nbasis,s,0,qq)
               write(iwr,'(a)') 
     &               ' checking orthogonality doubly-singly+virtuals'
               call ortche(vcopy,ndoc,vcopy(1,ndoc+1),nsoc+nvir,
     &                     nbasis,s,1,qq)
               write(iwr,'(a)')' checking orthogonality singly-virtuals'
               call ortche(vcopy(1,ndoc+1),nsoc,vcopy(1,ndoc+nsoc+1),
     &                     nvir,nbasis,s,1,qq)
               write(iwr,'(a)') 
     &                  ' ---------- end of current set ----------' 
            end if
c.....
c.....      put the vectors back
c.....
            call vecsca(v,vcopy(1,1          ),nbasis,idoc,ndoc)
            call vecsca(v,vcopy(1,ndoc+1     ),nbasis,isoc,nsoc)
            call vecsca(v,vcopy(1,ndoc+nsoc+1),nbasis,ivir,nvir)
c.....
220      continue
230   continue
240   continue
c.....
c.....reestablish correspondence of excitations per equivalence group
c.....this must be done for safety (signs etc.), things can get spoiled
c.....(numerically) because of the orthogonalisations. The check can
c.....terminate this run (crash = .true.)
c.....
      call reexc(hmoao,nao,nmo,v,iocvi,crash)
      if (crash) iprins = 100000
c
c..   locate possible zero norms, and virtuals not referred to.
c..   (e.g. pi-x pi-y for say Li2)
c..
      call zerono(v(1,nscf+1),s,nbasis,nzero,none,iredund)
c...  do not worry about zero virtuals which are zero because of harmon
      nzeroh = nzero
      if (oharm) nzeroh = nzeroh - (newbas1-newbas0)
      if (nzeroh.ne.0) then
         write(iwr,999) nzeroh,nzero
999      format(1x,i3,' virtuals have a (near)-zero norm of which ',
     1          i3,' due to harmon')
      end if
      if (none.gt.ncore.and.nitscf.le.1) then
         write(iwr,1111) none-ncore
1111     format(' # virtuals not referenced to :',i3,/)
      end if
c.....collect redundant virtuals, and sort in ascending order
c.....doubly and singly occupieds are added later on
c
      do i=1,ncore
         iredund(nzero+none+i) = iocvi(i)
      end do
c
c.....be sure not to count the core orbitals twice, for example:
c.....ethylene with in-plane orbitals as core, these are already
c.....counted as virtuals not referenced to.
c
      itwice=0
      do 2000 i=1,ncore
         j=iocvi(i)
         k=locati(j,iredund,nzero+none)
         if (k.ne.0) then
            itwice=itwice+1
            idref(itwice)=k
         end if
2000  continue
c
      if (itwice.gt.0) then
c
c...     get rid of double orbitals
c
         call ibubbl(idref,itwice,+1)
c        write (iwr,2122) itwice
c2122     format(' Number of redundant orbitals mentioned twice :'
c     &            ,i3,/)
         do k=1,itwice
            do i=idref(k),none+nzero+ncore-1
               if (i.gt.mxorbvb) call caserr('iredund overflow')
               iredund(i)=iredund(i+1)
            end do
c...        iredund has just been adjusted / adjust idref as well
            do i=k+1,itwice
               idref(i) = idref(i) - 1
            end do
            none = none - 1
         end do
      end if
      nredund = none + nzero + ncore
c
c...  remove redundancies due to doubly occupieds
c...  molecule like NO (2PHI)
c...  xn = equivalent with yn
c...  xo = equivalent with yo
c...  xn => xo is now removed (xo is doubly occupied)
c...  yn => yo is not removed (yo is singly occupied)
c...  if nothing is done, the symmetry will get lost.
c...  restore the excitation xn => xo to prevent this
c...  the orbital xo will get a minus sign in inforb.
c...  does not interfere with doubly/doubl2/redorb
c
      i=1
      do 2130 j=1,nequi
         if (iequi(j).eq.1) goto 2130
         do 2120 k=i,i+iequi(j)-2
            do 2110 l=k+1,i+iequi(j)-1
               do 2100 m=2,nscf+1
                  kk=locati(inteq(m,k),inforb(1,k),ninfor(k))
                  ll=locati(inteq(m,l),inforb(1,l),ninfor(l))
                  if (kk.ne.0.and.ll.eq.0) inforb(kk,k) = -inforb(kk,k)
                  if (ll.ne.0.and.kk.eq.0) inforb(ll,l) = -inforb(ll,l)
2100           continue
2110        continue
2120     continue
2130  i=i+iequi(j)
c
      do 260 i=1,nscf
         do 250 j=1,ninfor(i)
            i1  = iococ(iabs(inforb(j,i))+ncore,'loc')
            iorb = iocvi(i1)
            i2 = locati( iorb,iex(2,i),nex(i) )
            if (i2.ne.0) then
               i3 = locati(iorb,iredund,nredund)
               if (i3.eq.0) then
                  do k=1,nex(i)-i2
                     iex(i2+k,i) = iex(i2+k+1,i)
                     ieqsig(i2+k,i) = ieqsig(i2+k+1,i)
                  end do
                  nex(i) = nex(i) - 1
               end if
            end if
250      continue
         do 256 j=1,nredund
            i1 = locati( iredund(j),iex(2,i),nex(i) )
            if (i1.ne.0) then
               do k=1,nex(i)-i1
                  iex(i1+k,i) = iex(i1+k+1,i)
                  ieqsig(i1+k,i) = ieqsig(i1+k+1,i)
               end do
               nex(i) = nex(i) - 1
            end if
256      continue
         do 257 j=1,nforbi
            if (iex(1,i).eq.iforbi(j,1)) then
               ip = iococ(iforbi(j,2),'loc')
               if (ip.le.ncore) ip = 0
               if (ip.ne.0) then
                  ifo = iocvi(ip)
                  i1 = locati( ifo,iex(1,i),nex(i)+1)
                  if (i1.ne.0) then
                     do k=1,nex(i)-i1+1
                        iex(i1+k-1,i) = iex(i1+k,i)
                        ieqsig(i1+k-1,i) = ieqsig(i1+k,i)
                     end do
                     nex(i) = nex(i) - 1
                  end if
               end if
            end if
257      continue
260   continue
c
c...  singly occupieds are not oke yet (Koos/Anthony)
c...  doubly are not oke in case of NO example above.
c
      if (.not.super) then
         do 300 j=1,nscf
           do 290 i=2,nex(j)+1
               do 280 k=ncore        +1,ncore+ndoubly+nsingly
cKoos/Anthony  do 280 k=ncore+ndoubly+1,ncore+ndoubly+nsingly
                  if (iex(i,j).eq.iocvi(k)) 
     1                iex(i,j) = iococ(k,'ord') - ncore
280            continue
290         continue
300      continue
      end if
c
c...  now collect all redundant virtuals
c
      if (super) then
         nadd = ndoubly
      else
         nadd = ndoubly + nsingly
      end if
      do 1000 i=1,nadd
         iorb = iocvi(ncore+i)
         ithere = locati(iorb,iredund,nredund)
         if (ithere.eq.0) then
            nredund = nredund + 1
            iredund(nredund) = iorb
         end if
1000  continue
      call ibubbl(iredund,nredund,-1)
      if (iprins.ge.100) then
         write(iwr,333)ncore,ndoubly,nsingly,nactiv
333      format(' correspondence between internal orbitals and virtuals'
     &        ,' (after projections) core:',i3,' doubly:',i3,' singly:',
     &          i3,' activ:',i3)
         write(iwr,444) (iococ(i,'ord'),i=1,ncore+ndoubly+nsingly)
         write(iwr,444) (iocvi(i)-nactiv,i=1,ncore+ndoubly+nsingly)
444      format(1x,39i3)
      end if
c
c...  remove redundant excitations
c
      do 330 j=1,nscf
         do 320 i=2,nex(j)+1
            do 310 k=1,nredund
               if (iex(i,j).ge.iredund(k)) then
                  iex(i,j) = iex(i,j) - 1
               end if
310         continue
320      continue
330   continue

      nsa = nactiv + nbasis
      do 340 i=1,nredund
         call fmove(v(1,iredund(i)+1),v(1,iredund(i)),
     &                    (nsa-iredund(i)-i+1)*nbasis)
340   continue
      do 543 i=1,nredund
543   iredund(i) = iredund(i) - nscf
c
      if (iprins.gt.9999) then
         write(iwr,555) (iredund(ij),ij=1,nredund)
555      format(/' redundant atomic orbitals :',(t28,20i3))
         write(iwr,666)
666      format(' eq, stripped excitation patterns :')
         itel = 0
         ieq = 1
         do 350 i=1,nscf
            itel = itel + 1
            if (itel.gt.iequi(ieq)) then
               itel = 1
               ieq = ieq + 1
            end if
350      write(iwr,777) ieq,(iex(j,i)*ieqsig(j,i),j=1,nex(i)+1)
777      format(1x,i2,(t4,30i3))
         write(iwr,888)
888      format(' corresponding vectors :')
         call prvc(v,nsa-nredund,nbasis,v,'o','l')
      end if
      if (crash) call vberr('reexc failed because of lack of symmetry')
c
c...  if super cas remove active => active excitations
c...  if super act remove all excitations to virtuals
c....     (only here and in input referenced)
c
      if (super_cas.or.super_act) then
         do j=1,nscf
            ithere = locati(iex(1,j),isingly,nsingly)
700         do i=2,nex(j)+1
               jthere = locati(iex(i,j),isingly,nsingly)
               if ((jthere.ne.0.and.ithere.ne.0.and.super_cas)
     1         .or.(super_act.and.iex(i,j).gt.nscf)) then
                  nex(j) = nex(j) - 1
                  do k=i,nex(j)+1
                     iex(k,j)=iex(k+1,j)
                  end do
                  go to 700
               end if
            end do
         end do
      end if
c
      return
      end

c***********************************************************************
      function mvocvi(v,vec,nbasis,nvec,s,w)
c.....
c.....find virtual vector having the largest overlap with (occupied) v
c.....
      implicit REAL (a-h,o-z), integer(i-n)
INCLUDE(common/turtleparam)
INCLUDE(../m4/common/sizes)
INCLUDE(common/twice)
INCLUDE(common/splice)
INCLUDE(common/vbcri)
      dimension v(nbasis),vec(nbasis,nvec),s(*),w(nbasis)
c 
      smax = 0.0d0
      imax = 0
      call cntrc(s,v,w,nbasis,crilow)
      do 10 i=1,nvec
         overlap = dabs( ddot(nbasis,w,1,vec(1,i),1) )
         if (overlap.gt.smax) then
            smax = overlap
            imax = i
         end if
10    continue
      mvocvi = imax
      if (smax.lt.critor) then
         call caserr('maximum overlap is zero in mvocvi')
      endif
      return
      end
      integer function iococ(i,action)
      implicit none
      character(*) action
c.....
c..... get occupied order number now including ncore 
c..... this was iocvi(i,1)
c.....
INCLUDE(common/turtleparam)
INCLUDE(common/twice)
INCLUDE(common/splice)
      integer i,ii,k
      if (action.eq.'ord') then
         if (i.le.ncore) then
            iococ = i
         else if (i.gt.ncore.and.i.le.ncore+ndoubly) then
            iococ =idoubly(i-ncore) + ncore
         else
            iococ =isingly(i-ncore-ndoubly) + ncore
         end if
      else if (action.eq.'loc') then
         if (i.le.ncore) then
            iococ = i
            return
         else
            do k=1,ndoubly
               if (idoubly(k)+ncore.eq.i) then
                  iococ = k+ncore
                  return
               end if
            end do
            do k=1,nsingly
               if (isingly(k)+ncore.eq.i) then
                  iococ = k+ncore+ndoubly
                  return
               end if
            end do
         end if
         iococ = 0
         call caserr('no occupied orbitals found')
      else
         call caserr('illegal action in iococ')
      end if
      return
      end
 
c***********************************************************************
      subroutine ortche(a,na,b,nb,nbasis,s,ipar)
c
      implicit REAL  (a-h,o-z) , integer   (i-n)
c
c.....
c.....debugging routine. prints non-orthogonalities and traces
c.....overcomplete sets
c.....
INCLUDE(../m4/common/iofile)
INCLUDE(common/vbcri)
INCLUDE(../m4/common/gmempara)
INCLUDE(../m4/common/vcore) 
      dimension a(nbasis,na),b(nbasis,nb),s(*)
c
      kw  = igmem_alloc_inf(nbasis,'vbscf.m','ortche','kw',IGMEM_DEBUG)
      do 30 i=1,na
         call cntrc(s,a(1,i),Q(kw),nbasis,critor/100.0d0)
         if (ipar.eq.0) then
c....       a and b are the same set
            jstart = i + 1
         else
            jstart = 1
         end if
         do 20 j=jstart,nb
            if (ipar.eq.0) then
               aa = ddot(nbasis,Q(kw),1,a(1,j),1)
            else
               aa = ddot(nbasis,Q(kw),1,b(1,j),1)
            end if
            if (i.eq.j.and.ipar.eq.0) goto 20
            if(dabs(aa).gt.critor)write(iwr,11)i,j,aa
11          format(1x,i3,' overlaps with ',i3,' :',f15.11)
20       continue
30    continue
      call gmem_free_inf(kw,'vbscf.m','ortche','kw')
c
      kvc = igmem_alloc_inf(nbasis*na,'vbscf.m','ortche','kvc',
     @                      IGMEM_DEBUG)
      kw = igmem_alloc_inf(nbasis+na,'vbscf.m','ortche','kw',
     @                     IGMEM_DEBUG)
      iseta = idepen(a,Q(kvc),s,Q(kw),nbasis,na)
      write(iwr,22)
22    format(' first set ')
      call prvc(a,na,nbasis,a,'o','l')
      if (iseta.eq.0) then
         write(iwr,33)
33       format(' first set is linear independent')
      else
         write(iwr,44) iseta
44       format(' first set is linear dependent at vector ',i3)
      end if
      call gmem_free_inf(kw,'vbscf.m','ortche','kw')
      call gmem_free_inf(kvc,'vbscf.m','ortche','kvc')
      if (ipar.ne.0) then
         kvc = igmem_alloc_inf(nbasis*nb,'vbscf.m','ortche','kvc',
     @                         IGMEM_DEBUG)
         kw = igmem_alloc_inf(nbasis+nb,'vbscf.m','ortche','kw',
     @                        IGMEM_DEBUG)
         write(iwr,55)
55       format(' second set')
         call prvc(b,nb,nbasis,b,'o','l')
         isetb = idepen(b,Q(kvc),s,Q(kw),nbasis,nb)
         if (isetb.eq.0) then
            write(iwr,66)
66          format(' second set is linear independent')
         else
            write(iwr,77) isetb
77          format(' second set is linear dependent at vector ',i3)
         end if
         call gmem_free_inf(kw,'vbscf.m','ortche','kw')
         call gmem_free_inf(kvc,'vbscf.m','ortche','kvc')
      end if
c
      return
      end

c***********************************************************************
      integer function idepen(vcbu,vc,s,w,ndim,nmdim)
c
c     check on linear dependency by means of schmidt
c               return :
c                              0: no lin. dependency
c     any other positive integer: lin. dependence first noticed in this
c                                                              vector
c     main part stolen from normvc
c...  schmidt othogonalisation of nmdim vectors of length ndim
c
c..   vc  ndim*nmdim
c..   s   ndim*(ndim+1)/2
c..   w   ndim*2
c
      implicit REAL (a-h,o-z), integer (i-n)
      dimension vc(*),s(*),w(*)
INCLUDE(common/vbcri)
c
      call fmove(vcbu,vc,ndim*nmdim)
      idepen = 0
      if (nmdim.le.0) return
c
      if (ndim .eq. 1) then
         if (vc(1) * vc(1) .le. cridep) then
            idepen = 1
            return
         end if
         return
      end if
      np1 = ndim+1
      ivci = 1
      do 150 i=1,nmdim
         call cntrc(s,vc(ivci),w,ndim,critor/100.0d0)
         t = 0
         tnorm = 0
         ivcj = 1
         ivcw = ndim
         do 120 j=1,i
            tnorm = tnorm - t*t
            t = 0
            do 110 jj=1,ndim
               t = vc(ivcj)*w(jj) + t
  110       ivcj = ivcj + 1
            ivcw = ivcw + 1
  120    w(ivcw) = -t
         if (tnorm + t .le. cridep) then
            idepen = i
            return
         end if
         tnorm = 1/dsqrt(tnorm+t)
         w(ivcw) = 1
         do 140 k = 1,ndim
            t = 0
            ivcw = np1
            do 130 j=k,ivci,ndim
               t = vc(j)*w(ivcw) + t
  130       ivcw = ivcw + 1
            vc(ivci) = t*tnorm
  140    ivci = ivci + 1
  150 continue
c
      return
      end

c***********************************************************************
      subroutine proorb(v,s,nbasis,ndoc,nsoc)
c
      implicit REAL (a-h,o-z), integer (i-n)
c
      dimension v(nbasis,*),s(*)
INCLUDE(common/hsinfo)
INCLUDE(common/vbcri)
INCLUDE(../m4/common/gmempara) 
INCLUDE(../m4/common/vcore) 
c
c...  do all orthogonalisations allowed for 
c
c...  version 01-09-92 : lot of orthogonalisations/projections done in the
c...  old version were no longer necessary due to the special construction
c...  of the virtual space
c...  if (super) the virtuals are not orthogonal to (all of) the singles !
c
c...  lowdin orthogonalise the doubles
c
      if (ndoc.gt.1) then
         ks = igmem_alloc_inf(ndoc*(ndoc+1)/2,'vbscf.m','proorb',
     &                        'ks',IGMEM_DEBUG)
         kcmsq = igmem_alloc_inf(ndoc*ndoc,'vbscf.m','proorb',
     &                           'kcmsq',IGMEM_DEBUG)
         klambda = igmem_alloc_inf(ndoc,'vbscf.m','proorb',
     &                             'klambda',IGMEM_DEBUG)
         kvecmt = igmem_alloc_inf(nbasis,'vbscf.m','proorb',
     &                            'kvecmt',IGMEM_DEBUG)
         call lowdins(v,s,Q(ks),Q(kcmsq),
     &                Q(klambda),Q(kvecmt),
     &                ndoc,nbasis,crilow)
         call gmem_free_inf(kvecmt,'vbscf.m','proorb','kvecmt')
         call gmem_free_inf(klambda,'vbscf.m','proorb','klambda')
         call gmem_free_inf(kcmsq,'vbscf.m','proorb','kcmsq')
         call gmem_free_inf(ks,'vbscf.m','proorb','ks')
      end if
c
c...  schmidt orthogonalise the singles on the doubles
c
      kw = igmem_alloc_inf(nbasis,'vbscf.m','proorb',
     &                     'kw',IGMEM_DEBUG)
      call schmids(v(1,ndoc+1),v,s,Q(kw),nsoc,ndoc,
     &             nbasis,cridep,'norm')
      call gmem_free_inf(kw,'vbscf.m','proorb','kw')
c
      return
      end

c***********************************************************************
      subroutine reexc(haomo,nao,nmo,vec,iocvi,crash)
c
c.....reestablish the order and sign of (equivalent) excitations
c                                                                 
      implicit REAL (a-h,o-z), integer (i-n)
INCLUDE(common/turtleparam)
INCLUDE(../m4/common/sizes)
INCLUDE(common/infato)
INCLUDE(common/splice)
INCLUDE(common/brill)
      common /fast/ nospin,nosymm,nosymc
      logical nospin,nosymm,nosymc
INCLUDE(common/scftvb)
INCLUDE(common/vbcri)
INCLUDE(../m4/common/iofile)
INCLUDE(../m4/common/gmempara) 
INCLUDE(../m4/common/vcore) 
      logical internal
      dimension haomo(*),vec(nao,*),iocvi(*)
cjvl      dimension iclass(maxact)
      logical crash
      REAL lowesth1
      ind(i,j)=max(i,j)*(max(i,j)-1)/2+min(i,j)
c
c.....be sure to have a proper haomo
c...  vec is de active orbitals(nmo) and all orbtals (core+active+virtual)
c
      khao = igmem_alloc_inf(nao*(nao+1)/2,'vbscf.m','reexc',
     &                       'khao',IGMEM_DEBUG)
      khao2 = igmem_alloc_inf(nao*(nao+1)/2,'vbscf.m','reexc',
     &                       'khao2',IGMEM_DEBUG)
      call get1e(Q(khao),dummy,'v',Q(khao2))
      call gmem_free_inf(khao2,'vbscf.m','reexc','khao2')
      if (hybry) call clmetr(Q(khao),nao,'any')
      nam=nao+nmo
      kscr2 = igmem_alloc_inf(nam*nam+nam*nao,'vbscf.m','reexc',
     &                        'kscr2',IGMEM_DEBUG)
      kscr3 = igmem_alloc_inf(nam*nao,'vbscf.m','reexc',
     &                        'kscr3',IGMEM_DEBUG)
      call mult11(Q(khao),haomo,Q(kscr2),
     &            nam,nao,vec,Q(kscr3))
c
c...  as core orbitals may be anything 
c...  clear their interaction to avoid hybry trouble
c
      do i=1,ncore
         do j=1,nam
            haomo(ind(nmo+i,j)) = 0.0d0
         end do
      end do
c...  also clear interaction with yourself
c      do i=1,nmo
c         haomo(ind(i,i+nmo+ncore)) = 0.0d0
c      end do
c
c...  vec is de active orbitals(nmo) and all orbtals (core+active+virtual)
c...  check
c
c
      call gmem_free_inf(kscr3,'vbscf.m','reexc','kscr3')
      call gmem_free_inf(kscr2,'vbscf.m','reexc','kscr2')
      call gmem_free_inf(khao,'vbscf.m','reexc','khao')
c
      khaomo = igmem_alloc_inf(nam*(nam+1)/2,'vbscf.m','reexc',
     &                        'khaomo',IGMEM_DEBUG)
      call fmove(haomo,Q(khaomo),nam*(nam+1)/2)
      if (nosymm) goto 666
c
c.....next only works if equivalent orbitals are neighbours
c.....make an array that tells to which equivalence group a
c.....mo belongs
c
cjvl      i=0
cjvl      do 6 j=1,nequi
cjvl         do 5 k=1,iequi(j)
cjvl            i=i+1
cjvl            iclass(i)=j
cjvl5        continue
cjvl6     continue
c
      do 10 i=1,nmo
         iex(1,i)=i
10    continue
      lowesth1=9999999.9
      iorb1=1
      do 20 i=1,nequi
         iii=1
c...      apply hand given excitations if any
         if (nexeq(i).gt.0) then
            nn = nexeq(i)/(2*iequi(i))
            ii = ieqs(i)
            do j=1,nn
               iii = iii + 1
               do k=1,iequi(i)
                  iorbn = iexeq(ii)
                  im = iocvi(iexeq(ii+1)+ncore) 
                  iex(iii,iorbn) = im
                  Q(khaomo-1+ind(iorbn,im)) = 0.0d0
                  ii = ii + 2
               end do
            end do
c...           add occ-occ excitations
            iii = iii + 1
            do k=1,iequi(i)
               iorbn = iorb1 + k-1
               im = iocvi(iorbn+ncore) 
               iex(iii,iorbn) = im
               Q(khaomo-1+ind(iorbn,im)) = 0.0d0
            end do
         end if
c
         do 30 j=nmo+1,nmo+nao
            if (idequi.eq.2) then
               call abmax(Q(khaomo),
     &                    nmo+nao,nmo+1,'tri',iorb1,am,im)
            else
               im = j
            end if
            internal = im.le.2*nmo+ncore
            h1=dabs(Q(khaomo-1+ind(iorb1,im)))
            if (h1.gt.remcri) then
               if (h1.lt.lowesth1) lowesth1 = h1
               iii=iii+1
               iex(iii,iorb1)=im
               Q(khaomo-1+ind(iorb1,im)) = 0.0d0
               do 40 k = 2,iequi(i)
                  iorb2 = iorb1+k-1
                  if (idequi.eq.2) then
                    naomo = nmo + nao
                    if (im.le.2*nmo+ncore) naomo = 2*nmo+ncore
                    call abmax(Q(khaomo),naomo,nmo+1,'tri',iorb2,
     +                         bm,lm)
                  else
c                    lm = j-1
                     lm = nmo
22                   lm = lm + 1
*
*     what is the difference bewteen following two if statements?
*     the intel's idiot gets confused with the commented one. (zahid)
c                     if (dabs(Q(khaomo-1+ind(iorb2,lm))).lt.remcri
c     +                  .and.idequi.ne.3) goto 22
*
                     if (dabs(Q(khaomo-1+ind(iorb2,lm))).lt.
     +                  remcri.and.idequi.ne.3) goto 22
                  end if
                  h2=dabs(Q(khaomo-1+ind(iorb2,lm)))
cjvl                     if (locati(im,iocvi,nmo).ne.0) then
c.....             this is an internal excitation, be sure n
c.....             belongs to the same equivalence group
c.....             (had problems with nh3 symmetry)
c.....             how does this behave with super option?
cjvl               disabled ; don't understajnd class
cjvl                     iom = locati(im,iocvi,nmo)-ncore
cjvl                     iol = locati(lm,iocvi,nmo)-ncore
cjvl                     if (iclass(iom).ne.iclass(iol)) then
cjvl                        write(iwr,*)im,lm,iom,iol,nmo,ncore,k,i,iequi(i)
cjvl                        write(iwr,*) 'classes ',
cjvl     1                              (iclass(l),l=1,max(im,lm)+nmo+ncore)
cjvl                        call caserr('class ?')
cjvl                     end if
cjvl                  end if
                  iex(iii,iorb2)=lm
                  Q(khaomo-1+ind(iorb2,lm)) = 0.0d0
40             continue
            end if
30       continue
         if ((iii-1).ne.nex(iorb1)) then
            write(iwr,*) 'Number of excitations has changed for',iorb1
            write(iwr,621) nex(iorb1),iii-1
            do jj = 1, iequi(i)
               write(iwr,*) (iex(kk,iorb1+jj-1),kk=1,iii)
            end do
621         format(/,' number before: ',i3,'. number after: ',i3)
            write(iwr,622) lowesth1
622         format(/,' Lowest H-matrix value: ',e14.7)
            write(iwr,*) 'Solve this by either switching off symmetry'
            write(iwr,*) 'or by using "remove" directive in scf'
            call caserr('Number of excitations has changed')
         end if
         iorb1=iorb1+iequi(i)
20    continue
c
666   call gmem_free_inf(khaomo,'vbscf.m','reexc','khaomo')
      i=0
      do 170 j=1,nequi
         do 160 k=1,iequi(j)
            i=i+1
            do 150 l=2,nex(i)+1
               h1=haomo(ind(iex(1,i),iex(l,i)))
               ieqsig(l,i)=-int(dsign(1.0d0,h1))
150         continue
160      continue
170   continue
c
      return
      end

c***********************************************************************
      subroutine subex(iorb1,idoc,ndoc,isoc,nsoc,ivir,nvir,iex2,nbasis)
      implicit REAL (a-h,o-z), integer(i-n)
c.....
c.....find subsets in excitations
c.....
      dimension idoc(*),isoc(*),ivir(*),iex2(nbasis+1,nscf)
INCLUDE(common/turtleparam)
INCLUDE(../m4/common/sizes)
INCLUDE(common/twice)
INCLUDE(common/brill)
      logical there1,there2
c.....
c.....loop over active occupied orbitals :
c.....
      ndoc = 0
      nsoc = 0
      nvir = 0
      ithere = locati(iorb1,isingly,nsingly)
      if (ithere.ne.0) then
         nsoc = nsoc + 1
         isoc(nsoc) = iorb1
      else
         nvir = nvir + 1
         ivir(nvir) = iorb1
      end if
c.....
      do 80 i=1,nscf
c.....
c.....   loop over sets
c.....
         do 50 k=2,nexw(i)+1
            iorb2 = iex2(k,i)
            if (iorb2.ne.0.and.iorb1.ne.iorb2) then
               do 20 l=i,nscf
                  there1 = .false.
                  there2 = .false.
                  do 10 m=2,nexw(l)+1
                     if (iex2(m,l).eq.iorb1) there1 = .true.
                     if (iex2(m,l).eq.iorb2) there2 = .true.
10                continue
                  if ((there1.and..not.there2).or.
     &                (there2.and..not.there1)   ) goto 50
20             continue
c.....
c.....         apparently orb1 and orb2 always appear together
c.....
               ithere = locati(iorb2,isingly,nsingly)
               if (ithere.ne.0) then
                  nsoc = nsoc + 1
                  isoc(nsoc) = iorb2
               else
                  nvir = nvir + 1
                  ivir(nvir) = iorb2
               end if
c.....
c.....         clear redundant information
c.....
               do 40 m=1,nscf
                  do 30 n=2,nexw(m)+1
                     if (iex2(n,m).eq.iorb2) iex2(n,m) = 0
30                continue
40             continue
            end if
50       continue
80    continue
c.....
c.....find occupieds, classify them and remove iorb1
c.....
      do 100 i=1,nscf
         do 90 j=2,nexw(i) + 1
            if (iex2(j,i).eq.iorb1) then
               ithere = locati(iex2(1,i),idoubly,ndoubly)
               if (ithere.ne.0) then
                  ndoc = ndoc + 1
                  idoc(ndoc) = iex2(1,i)
               end if
               ithere = locati(iex2(1,i),isingly,nsingly)
               if (ithere.ne.0) then
                  iorb = isingly(ithere)
                  ither2 = locati(iorb,isoc,nsoc)
                  if (ither2.eq.0) then
                     nsoc = nsoc + 1
                     isoc(nsoc) = iorb
                  end if
               end if
               iex2(j,i) = 0
            end if
90       continue
100   continue
      return
      end

c***********************************************************************
      subroutine zerono(v,s,nbasis,nzero,none,izero)
      implicit REAL (a-h,o-z), integer(i-n)
c
c..   find out about zero norms, and redundant orbitals
c..   (i.e. not refered to in excitations)
c
INCLUDE(common/turtleparam)
INCLUDE(../m4/common/sizes)
INCLUDE(common/brill)
INCLUDE(common/vbcri)
      dimension v(nbasis,nbasis),s(*),izero(*)
      nzero = 0
      none = 0
      do 20 i=1,nbasis
         vnor = vnormvb(v(1,i),s,nbasis)
         if (vnor.lt.crinor) then
            nzero = nzero + 1
            izero(nzero+none) = i + nscf
            goto 20
         end if
         do 10 j=1,nscf
            ithere = locati(i+nscf,iex(2,j),nex(j))
            if (ithere.ne.0) goto 20
10       continue
         none = none + 1
         izero(nzero+none) = i + nscf
20    continue
      return
      end

c***********************************************************************
      function vnormvb(v,s,nbasis)
      implicit REAL (a-h,o-z), integer(i-n)
      dimension v(nbasis),s(*)
      vnormvb = 0.0d0
      do 30 i=1,nbasis
         ii = i*(i-1)/2
         vv = v(i)
         do 10 j=1,i
10       vnormvb = vnormvb + vv * s(ii+j) * v(j)
         do 20 j=i+1,nbasis
20       vnormvb = vnormvb + vv * s(j*(j-1)/2+i) * v(j)
30    continue
      return
      end
 
c***********************************************************************
      subroutine reorb(idet,ndets,nelec,nalfa,iorth,iscr,coeff,
     &                                          idetps,ndettot)
c
      implicit REAL  (a-h,o-z) , integer   (i-n)
c
c.....
c.....in this routine the determinants in idet are reordered so that
c.....the orbitals appear in mutually orthogonal groups, for the benifit
c.....of the matrix element evaluation. this is done for alpha and beta
c.....spin separately, as spin introduces an extra orthogonality on top
c.....of the spatial one. in the i-th position iorth must contain a
c.....number between 1 and the number of orthogonality groups, that is
c.....characteristic for the i-th orbital its group. parity changes are
c.....accounted for via the coefficients of the determinants.
c.....
INCLUDE(common/turtleparam)
      dimension idet(ndets*nelec),iorth(*),iscr(*),coeff(*),idetps(*)
      istart = 1
      it     = 0
      do 50 m=1,ndets
         ipar  = 1
         i     = istart
         irrep = 1
10       do 20 k=i,istart+nalfa-1
            if (irrep.eq.iorth(idet(k))) then
               iii      =  idet(k)
               idet(k)  =  idet(i)
               idet(i)  =  iii
               if (i.ne.k) ipar = -ipar
               i        =  i + 1
            end if
20       continue
         if (i-istart+1.le.nalfa) then
            irrep = irrep + 1
            if (irrep.gt.mxorbvb) then
               call caserr('problem in reorb')
            else
               goto 10
            endif
         end if
         irrep  = 1
         istart = istart + nalfa
         i      = istart
30       do 40 k=i,istart+nelec-nalfa-1
            if (irrep.eq.iorth(idet(k))) then
               iii     =  idet(k)
               idet(k) =  idet(i)
               idet(i) =  iii
               if(i.ne.k) ipar = -ipar
               i       =  i + 1
            end if
40       continue
         if (i-istart+1.le.nelec-nalfa) then
            irrep = irrep + 1
            if (irrep.gt.mxorbvb) then
               call caserr('problem in reorb')
            else
               goto 30
            endif
         end if
         istart = istart + nelec - nalfa
         if(ipar.eq.-1) then
            it = it + 1
            iscr(it) = m
         end if
50    continue
      do 70 i=1,it
         do 60 j=1,ndettot
            if(idetps(j).eq.iscr(i)) coeff(j) = -coeff(j)
60       continue
70    continue
      return
      end

c***********************************************************************
      subroutine schmids(v1, v2, s, w, n1, n2, n, critlow,norm)
c
      implicit REAL  (a-h,o-z) , integer   (i-n)
c
c
c     schmidt orthogonalization of n1 vectors(v1) of dimension n on
c     n2 vectors(v2) of dimension n
c     the resulting vector-set is not internally orthogonalized
c     results are returned in v1
c
c     version with metric matrix
c
      dimension v1(*),v2(*),s(*),w(*)
      character*(*) norm
c
      if (n1.le.0.or.n2.le.0) return
c
      k1 = -n
      do 60 i1=1,n1
        k1 = k1 + n
         call cntrc(s, v1(k1 + 1), w, n,critlow/1000.0d0)
         k2 = -n
         do 60 i2=1,n2
            k2 = k2 + n
            res = 0.0d0
            do 50 j=1,n
   50       res = res + w(j)*v2(k2+j)
            do 60 j=1,n
   60    v1(k1+j) = v1(k1+j) - res*v2(k2+j)
c
      if (norm.ne.'norm') return
c
c...  normalize set v1
c
      k1 = 0
      do 100 i1=1,n1
         call cntrc(s,v1(k1+1),w,n,critlow/1000.0d0)
         ss = ddot(n,v1(k1+1),1,w,1)
         if (ss.lt.critlow) then
c....     vector is effectively 0.0 / but who cares in vb
            tnorm = 0.0d0
         else
            tnorm = 1.0d0/dsqrt(ss)
         end if
         do 90 j=1,n
90       v1(k1+j) = v1(k1+j)*tnorm
100   k1 = k1 + n
c
      return
      end

c***********************************************************************
      subroutine schmidt(s,t,n,itroub)
c
      implicit REAL (a-h,o-z), integer (i-n)
c
      dimension s(*),t(n,n)
INCLUDE(common/vbcri)
      ind(i,j) = max(i,j) * (max(i,j)-1) / 2 + min(i,j)
      itroub = 0
      call vclr(t,1,n*n)
      do 50 i=1,n
         t(i,i) = 1.0d0
         if (s(ind(i,i)).lt.cridep) then
            itroub = i
            return
         end if
         do 30 j=i+1,n
            t(i,j) = -s(ind(i,j))/s(ind(i,i))
            s(ind(j,j)) = s(ind(j,j)) + t(i,j) * s(ind(i,j))
            do 10 k=1,i-1
               t(k,j) = t(k,j) + t(k,i)*t(i,j)
10          continue
            do 20 k=j+1,n
               s(ind(k,j)) = s(ind(k,j)) + s(ind(k,i))*t(i,j)
20          continue
30       continue
         root = 1/dsqrt(s(ind(i,i)))
         do 40 j=1,i
            t(j,i) = t(j,i) * root
40       continue
50    continue
      return
      end

c***********************************************************************
      subroutine sigorb(s,n)
c
      implicit REAL (a-h,o-z), integer (i-n)
c
c...  prints the sign matrix of triangle s
c
      dimension s(*),iscr(100),scr(100)
      character*8 scr
INCLUDE(common/vbcri)
INCLUDE(../m4/common/iofile)
      ind(i,j) = max(i,j) * (max(i,j)-1)/2 + min(i,j)
c
      if (n.gt.100) then
         write(iwr,1) 
1        format(' ** SIGORB cannot print; adjust dimensions')
         return
      end if
c
      do 10 i=1,n
10    iscr(i) = i
c
      write(iwr,20) (iscr(ij),ij=1,n)
20    format(/,' sign matrix between the ao''s: ',/,
     &                                   (t4,40i3),/)
c
      do 50 i=1,n
         do 30 j=1,n
            if (s(ind(i,j)).gt.critor) then
               scr(j) = '+'
            else if (s(ind(i,j)).lt.-critor) then
               scr(j) = '-'
            else
               scr(j) = ' '
            end if
30       continue
         write(iwr,40)i,(scr(ij),ij=1,n)
40       format(1x,i3,1x,(40a3))
50    continue
c
      return
      end

c***********************************************************************
      subroutine sym1(hmat,ndim,isy,nrep,mode)
c
      implicit none
c
c...  find symmetry classification from 1-electron matrix
c...  results (symmetry numbers) are returned in isy
c
INCLUDE(common/turtleparam)
INCLUDE(../m4/common/sizes)
INCLUDE(../m4/common/iofile)
INCLUDE(common/vbcri)
      REAL hmat(*)
      integer nrep,ndim,isy(ndim)
      integer ig(mxorbvb * 10)
      character *(*) mode
      integer i,iky,irf,irt,k,n,ij,j
c
      iky(i) = i*(i-1)/2
c
      call setsto(ndim,0,isy)
      nrep = 0
c
      do 30 i=1,ndim
         if (isy(i).eq.0) then
            nrep = nrep + 1
            isy(i) = nrep
         end if
         do 20 j=i+1,ndim
            if (dabs(hmat(iky(j)+i)).gt.critor) then
               if (isy(j).gt.0.and.isy(j).ne.isy(i)) then
                  irf = max(isy(i),isy(j))
                  irt = min(isy(i),isy(j))
                  nrep = nrep - 1
                  do 10 k=1,ndim
                     if (isy(k).eq.irf) isy(k) = irt
10                continue
               else
                  isy(j) = isy(i)
               end if
            end if
20       continue
30    continue
c
c.... print orthogonality information if requested
c
      if (mode.eq.'print') then
         if (nrep.eq.1) then
            write(iwr,77)
77          format(//' ',' all orbitals are within the same '
     &                ,'orthogonality class')
         else
            write(iwr,33) nrep
33       format(//' ',i4,' orthogonality classes have been recognised ',
     *          'in the current basisset')
            if (ndim.eq.nrep) then
               write(iwr,44)
44       format(/' ','this seems to be a completely orthogonal problem')
            else
               do 70 i=1,nrep
                  n = 0
                  do 60 j=1,ndim
                     if (isy(j).eq.i) then
                        n = n + 1
                        ig(n) = j
                     end if
60                continue
                  write(iwr,55) i,n,(ig(ij),ij=1,n)
55             format(/' ','group ',i3,' are ',i3,' orbitals :',50i3)
70             continue
            end if
         end if
      end if
      return
      end

c***********************************************************************
      subroutine toetoe(vec,haomo,nao,nmo)
c
c.....establish relevant excitations for the scf process by looking
c.....at one-electron integrals and find out about symmetry
c..... vec contains active orbitals followed by the rest (incl. active/core)
c
      implicit REAL (a-h,o-z),integer (i-n)
INCLUDE(common/turtleparam)
INCLUDE(../m4/common/sizes)
INCLUDE(../m4/common/vcore) 
INCLUDE(../m4/common/gmempara)
INCLUDE(../m4/common/gmemdata)
INCLUDE(common/splice)
INCLUDE(common/brill)
INCLUDE(common/scftvb)
INCLUDE(common/vbcri)
INCLUDE(common/infato)
INCLUDE(common/twice)
      common /fast/ nospin,nosymm,nosymc
      logical nospin,nosymm,nosymc
INCLUDE(../m4/common/iofile)
INCLUDE(common/vbproper)
      dimension vec(nao,*),haomo(*)
      logical symcha,eqfa,ochatao,oflop
      ind(i,j)=max(i,j)*(max(i,j)-1)/2+min(i,j)
c
c.....get one-electron integrals, remove kinetic component
c
      khao = igmem_alloc_inf(nao*(nao+1)/2,'vbscf.m','toetoe',
     &                       'khao',IGMEM_DEBUG)
      khao2 = igmem_alloc_inf(nao*(nao+1)/2,'vbscf.m','toetoe',
     &                       'khao2',IGMEM_DEBUG)
      call get1e(Q(khao),dummy,'v',Q(khao2))
      call gmem_free_inf(khao2,'vbscf.m','toetoe','khao2')
      if (hybry) call clmetr(Q(khao),nao,'any')
c
c.....construct matrix representation of the electron-nuclei
c.....attraction operator in basis of occupieds and virtuals
c
      nam=nao+nmo
      kscr2 = igmem_alloc_inf(nam*nam+nam*nao,'vbscf.m','toetoe',
     &                        'kscr2',IGMEM_DEBUG)
      kscr3 = igmem_alloc_inf(nam*nao,'vbscf.m','toetoe',
     &                        'kscr3',IGMEM_DEBUG)
      call mult11(Q(khao),haomo,Q(kscr2),nam,nao,vec,
     &            Q(kscr3))
c
c...  as core orbitals may be anything 
c...  clear their interaction to avoid hybry trouble
c
      do i=1,ncore
         do j=1,nam
            haomo(ind(nmo+i,j)) = 0.0d0
         end do
      end do
c
      call gmem_free_inf(kscr3,'vbscf.m','toetoe','kscr3')
      call gmem_free_inf(kscr2,'vbscf.m','toetoe','kscr2')
      call gmem_free_inf(khao,'vbscf.m','toetoe','khao')
c
c.....find out which scf orbitals are equivalent
c..... if we forced it bypass it
c
      if (idequi.ne.0) go to 51
      if (nosymm) then
c      if (nosymm.and.nequi.eq.0) then
         nequi=nscf
         do 10 i=1,nscf
            iequi(i)=1
10       continue
      else
         i=1
         nequi=0
         symcha=.false.
20       ii=ind(i,i)
         nequi=nequi+1
         iold=iequi(nequi)
         iequi(nequi)=1
         do 50 j=i+1,nscf
            jj=ind(j,j)
            if (dabs(1.0d0-haomo(jj)/haomo(ii)).lt.symcri) then
c.....         do second check (hopefully sufficient)
               if (iequi(nequi).eq.1) then
c.....            sumi not known yet
                  sumi=0.0d0
                  do 30 k=nmo+1,nam
                     ki=ind(k,i)
                     sumi=sumi+dabs(haomo(ki))
30                continue
               end if
               sumj=0.0d0
               do 40 k=nmo+1,nam
                  kj=ind(k,j)
                  sumj=sumj+dabs(haomo(kj))
40             continue
               if (dabs(1.0d0-sumj/sumi).lt.symcri) then
c.....            if i and j are equally occupied equivalency is "exact"
c.....            if not, check if equivalency is "forced"
                  id=locati(i,idoubly,ndoubly)
                  jd=locati(j,idoubly,ndoubly)
                  if (id.eq.0.neqv.jd.eq.0) then
c.....               occupancy i .ne. occupancy j
c.....               if equivalency is not imposed goto next j
                     if ((exact(nequi)).or.(.not.exact(nequi).and.
     &                   iold.eq.1)) goto 50
                  end if
c.....            i and j are found/forced to be equivalent
c.....            equivalent orbitals ougth to be neighbours
                  if (i+iequi(nequi).ne.j) then
                    write(iwr,501) i,j
501    format(' orbital',I4,' and',I4,' found to be equivalent')
                    call vberr('reordering required due to symmetry')
                  end if
c.....            would be nice if reordering is done automaticly
c.....            in the future (equiv then needs some changes!)
                  iequi(nequi)=iequi(nequi)+1
               end if
            end if
50       continue
c
c.....  koos: don't use symmetry in case an orbital appears alone (why?)
         if (iequi(nequi).eq.1) exact(nequi)=.false.
c
         if (iold.ne.0.and.iold.ne.iequi(nequi)) symcha=.true.
         i=i+iequi(nequi)
         if (i.le.nscf) goto 20
         if (symcha) then
            write(iwr,*) ' --> warning <-- symmetry has changed,'
            if (nequi.eq.nscf) then
               write(iwr,*) 
     &               '                 there are no equivalent',
     &               ' orbitals anymore, symmetry switched off'
               nosymm=.true.
            else
               do 55 i=1,nequi
                  write(iwr,551) i,iequi(i)
551   format ('                 the',I4,'-th equivalence group has',
     &        I4,' members now')
55             continue
            end if
         end if
      end if
51    continue
c
      if (super_hybrid) return
c
c.....find out which virtuals have the right symmetry to mix
c
      call izero(maxact*maxex,iexw,1)
      do 80 i=1,nscf
         if (ofrzmo(i)) go to 80
         k=0
         iexw(1,i)=i
c...      next part may be done in reexc but is not needed
         if (hybry) then
c.....      determine the atom orbital i belongs to
            iatom=0
            do 60 j=1,natom
               ithere=locati(i,iacat(1,j),nacat(j))
               if (ithere.ne.0) iatom=j
60          continue
            if (iatom.eq.0) call vberr('inconsistent atom definition')
         end if
         do 70 j=nmo+1,nmo+nao
            if (hybry) then
c.....         determine the atom orbital j belongs to
c.....          find largest coefficient to determine atom
               call abmax(vec(1,j),nao,1,'sq',1,am,im)
               if (.not.ochatao(iatom,im,'atomao')) goto 70
            end if
            if (dabs(haomo(ind(i,j))).gt.remcri
     &          .or.orunprop.or.onosym) then
               k=k+1
               iexw(k+1,i)=j
            end if
70       continue
         nexw(i)=k
80    continue
c
c...  If for a group of equivalent orbitals the equivalence is
c...  not exact make the equivalence of other groups false too.
c...  This is done for safety (see NO example in redorb: if on-
c...  ly xo & yo are forced the answer will be wrong).
c
      eqfa=.false.
      do 100 i=1,nequi
         if (iequi(i).gt.1) then
            if (.not.exact(i)) eqfa=.true.
         end if
100   continue
      if (eqfa) then
         do 110 i=1,nequi
            exact(i)=.false.
110      continue
      end if
c
      return
      end

c***********************************************************************
      subroutine vbper(ciint,bril,grh,grs,diagh,diags,igr_brill)
c
      implicit REAL  (a-h,o-z) , integer   (i-n)
c
c...  merge the brill-vector with perturbation estimates
c...  Bij = Hij0-E0*Sij0/(E0*Sijij - Hijij)
c
INCLUDE(common/scftvb)
INCLUDE(common/hsinfo)
INCLUDE(common/turtleparam)
INCLUDE(../m4/common/sizes)
INCLUDE(../m4/common/iofile)
INCLUDE(common/brill)
INCLUDE(common/vbpert)
c
      dimension ciint(*),bril(*),grh(*),grs(*),diagh(*),diags(*)
      integer igr_brill(5,*)
c.....
c.....get grh (first column of brill.matrix), grs (corresponding
c.....overlap), diagh (diagonal of brill) and diags (diagonal of s)
c.....
      call search(ipbloc,ihfile)
      call reads(grh,nbrill,ihfile)
      call reads(diagh,nbrill,ihfile)
      call reads(grs,nbrill,ihfile)
      call reads(diags,nbrill,ihfile)
c
c...  read excitation types from ed7 (in igr_brill)
c
      k7igr_brill = kscra7vb('k7igr_brill',nbrill*5,'i','r')
      call readi(igr_brill,nbrill*5,k7igr_brill,num8)
c
c...   check phase of BI-vector (PT assumes +)
c
      phase = 1.0d0*dsign(1.0d0,ciint(1))
      e0 = diagh(1) / diags(1) - shiftp
      bril(1) = ciint(1)
c
c...  grh = hijij - e0*sij0
c
      call daxpy(nstruc-1,-e0,grs(2),1,grh(2),1)
c
c...  diagh = eij - e0*sijij
c
      call daxpy(nstruc-1,-e0,diags(2),1,diagh(2),1)
c
c...  bij = -grh / diagh
c
      do 10 i=2,nstruc
         bril(i) = -phase * grh(i)/diagh(i)
10    continue
c
c...  now insert the variational values from ciint
c
      iv = 1
      ebi = 0.0d0
      do 20 ib=2,nstruc
         if (igr_brill(4,ib).eq.1) then
            iv = iv + 1
            bril(ib) = ciint(iv)
         else
            ebi = ebi + bril(ib) * grh(ib) * phase
         end if
20    continue
      return
      end

c***********************************************************************
      subroutine vecgat(vold,vnew,nbasis,map,new)
      implicit REAL (a-h,o-z), integer(i-n)
c.....
c.....gather vectors from vold in vnew according to map
c.....
      dimension vold(nbasis,*),vnew(nbasis,*),map(new)
      do 10 i=1,new
10    call fmove(vold(1,map(i)),vnew(1,i),nbasis)
      return
      end

c***********************************************************************
      subroutine vecsca(vnew,vold,nbasis,map,new)
      implicit REAL (a-h,o-z), integer (i-n)
      dimension vnew(nbasis,*),vold(nbasis,*),map(new)
      do 10 i=1,new
10    call fmove(vold(1,i),vnew(1,map(i)),nbasis)
      return
      end

c***********************************************************************
      subroutine virtual(occ,virt,ndim,ncore,nocc,nvirt,maxvirt)
c
      implicit REAL (a-h,o-z), integer(i-n)
c
c...  the virtual space is constructed (see thesis Verbeek, page 76)
c...  occup : occupied orbitals (input)
c...  virt  : virtual orbitals  (output)
c...  ndim = nbasis
c...  occ  and virt are consecutive parts of v
c...  both occupied and virtual orbitals may be canonicalised or localised
c...   in common/twice/ 
c...  ndoubly,idoubly : # always doubly occupied orbitals ; doc
c...  nsingly,isingly : # always doubly occupied orbitals ; voc
c
c...  in the virt array the occupied and virtuals are consecutive
c...  this was done so in the original virtual and weeded out by redorb
c...  we now  emulate this
c
c...  for dynamic memory allocation debugging purpose, contains IGMEM_ vars.
INCLUDE(../m4/common/gmemdata)
INCLUDE(../m4/common/gmempara)
INCLUDE(../m4/common/vcore)
c
INCLUDE(common/turtleparam)
INCLUDE(../m4/common/sizes)
INCLUDE(../m4/common/iofile)
INCLUDE(../m4/common/infoa)
INCLUDE(common/infato)
INCLUDE(common/vbequiv)
INCLUDE(common/scftvb)
INCLUDE(common/twice)
INCLUDE(common/vbvirt)
INCLUDE(common/vbcri)
INCLUDE(common/brill)
INCLUDE(common/tractlt)
INCLUDE(common/vbpert)
INCLUDE(../m4/common/harmon)
INCLUDE(../m4/common/mapper)
INCLUDE(../m4/common/tran)
c
      dimension occ(ndim,nocc),virt(ndim,maxvirt)
c
      logical frozen,normal
      character*5 harmvb
c
c...  now we can diagonalise in harmonic basis if required
c
c     itri(i,j) = max(i,j)*(max(i,j)-1)/2 + min(i,j)
c...  optional switches: (now set in vbin)
c...  idempotent = make a proper, idempotent projection operator
c...  canonicalise = make the virtuals reflect the point group symmetry of the
c...                 molecule
c...  aos = use atomic orbitals to mix (not tested thoroughly)
c...  super_hybrid : use exactly the hybrid excitations and 
c...                 thus those ao's as virtuals (in /twice/)
c...                 assume last one is most identical to mo and thus skipped
c
c     idempotent=.false.
c     canonicalise=1
c     aos=.false.
c...    frozen signals frozen atoms 
       do i=1,nbasis
         ilifq(i) = (i-1)*nbasis
       end do
c      do i=1,nbasis
c      iky(i)=i*(i-1)/2
c      end do
      if (ndim.ne.nbasis) call caserr('nbasis clash in virtual')
      if (nscf.ne.nocc-ncore) call caserr('nscf clash in virtual')
      if (nsingly.ne.nocc-ndoubly-ncore) 
     1   call caserr('nsingly clash in virtual')
c
      frozen = .false.
      do i=1,natom
        if (ifrzat(i).gt.1) frozen = .true.
      end do
c
      normal = .true.
      if (aos.or.super_hybrid) normal = .false.
c

      if (normal) then
c
c...    normal vbscf
c
_IF(debug)
         write(iwr,*) '((0)) VIRTUAL - AOS_SUPERHYBRID = FALSE __ OK'
_ENDIF
c
c...  get metric
c

         ksao = igmem_alloc_inf(nbasis*(nbasis+1)/2,'vbscf.m','virtual',
     &                          'ksao',IGMEM_DEBUG)
         kscr = igmem_alloc_inf(nbasis*(nbasis+1)/2,'vbscf.m','virtual',
     &                          'kscr',IGMEM_DEBUG)
         call get1e(Q(ksao),dummy,'s',Q(kscr))
         call gmem_free_inf(kscr,'vbscf.m','virtual','kscr')
         if (hybry) call clmetr(Q(ksao),nbasis,'any')
c
         koccp = igmem_alloc_inf(ndim*nocc,'vbscf.m','virtual',
     &                           'koccp',IGMEM_DEBUG)
         call fmove(occ,Q(koccp),nbasis*nocc)
c
c...  be sure to clear parts that were not cleared because of frozen
c...  or core ---
c
         if (clean.and.hybry) then
_IF(debug)
            write(iwr,*) '((1)) VIRTUAL - CLEAN/HYBRY __ OK'
_ENDIF
            call clvecvb(Q(koccp),
     &                   nocc,nbasis,'all','any','in virtual')
         end if
c
c...  make projection operator
c
         if (idempotent) then
_IF(debug)
            write(iwr,*) '((2)) VIRTUAL - INDEMPOTENT __ OK'
_ENDIF
c...  schmidt orthogonalise the occupieds
            kscr    = igmem_alloc_inf(ndim*2,'vbscf.m','virtual','kscr',
     &                                IGMEM_DEBUG)
            call normvc(Q(koccp),Q(ksao),
     &                  Q(kscr),nbasis,nocc,cridep)
            call gmem_free_inf(kscr,'vbscf.m','virtual','kscr')
         end if
c
c...  start of possible *harmonic*
c...  in vbfirst start orbitals were projected to harmonic
c...  at the end of this routine the zero vectors are move to the back
c...  elsewhere the check on zero vectors is suitably relaxed
c...  all is (of course) marked by oharm
c
         ns=newbas0*(newbas0+1)/2
         if (oharm) then
_IF(debug)
            write(iwr,*) '((3)) VIRTUAL - OHARM __ OK'
_ENDIF
            otran = .false.
            write(harmvb,'(a2,i3)') 'vb',nocc
            kscr = igmem_alloc_inf(nbasis*nbasis,'vbscf.m','virtual',
     &                             'kscr',IGMEM_DEBUG)
c
c...       get occupied orbitals in adapted basis
c
            call anorm(Q(kscr),Q(1))
            call tback(Q(koccp),ilifq,Q(koccp),ilifq,nocc)
c...     get vectors closely packed but in the 2-d-array
            call comharm(Q(koccp),harmvb,ilifq)
c
c...   get s-matrix in reduced harmonic basis
c
            call comharm(flop,'ctrans',flop)
            call tranp(Q(ksao),Q(kscr))
            call dcopy(ns,Q(kscr),1,Q(ksao),1)
            call gmem_free_inf(kscr,'vbscf.m','virtual','kscr')
         else
            if (newbas0.ne.nbasis) 
     1         call caserr('harmonic blows in virtual')
         end if
c
_IF(debug)
         write(iwr,*) '((_)) VIRTUAL - NO IF BLOCK AFTER (3) __ OK'
_ENDIF
c
c...    get the occ density matrix
c
         kprma = igmem_alloc_inf(ns,'vbscf.m','virtual','kprma',
     &                           IGMEM_DEBUG)
c
         call vclr(Q(kprma),1,ns)
         k=0
         do i=1,newbas0
            do j=1,i
               k = k + 1
               do n=1,nocc
                  if (super.and.n.gt.ncore) then
c...  exclude the singly occupieds, as they may be very similar
                     if (locati(n-ncore,isingly,nsingly).ne.0) exit
                  end if
c...  since comharm may have compacted vectors (taking ilifq along)
c          print *,'remco',j,i,n,nocc,'kprma',kprma+k-1,koccp-1,ilifq(n)
c            print *,'remcolocm',loc( Q(kprma+k-1)),
c.     &       loc(Q(koccp-1+i+ilifq(n))),
c     &  loc(Q(koccp-1+j+ilifq(n)))
c               call remcoi('prma again',k) 
                  Q(kprma+k-1) = Q(kprma+k-1) +
     &                  Q(koccp-1+i+ilifq(n))
     &                 *Q(koccp-1+j+ilifq(n))
               end do
            end do
         end do
c...  clear density matrix to conform to hybrid
c...  (*to be sure*  is *not* needed* ??? and confuses harmon
ccc      call clmetr(Q(kprma),nbasis,'small')
c
         if (iprins.ge.50000) then
            write(iwr,*) ' projection operator'
            call tripri(Q(kprma),newbas0)
         end if
c
c...  diagonalise the projection operator, virtuals have eigenvalue zero
c
         ns = nbasis*(nbasis+1)/2
         kvec  = igmem_alloc_inf(nbasis*nbasis,'vbscf.m','virtual',
     &                           'kvec',IGMEM_DEBUG)
         kval  = igmem_alloc_inf(nbasis,'vbscf.m','virtual','kval',
     &                           IGMEM_DEBUG)
         kscr    = igmem_alloc_inf(2*ns+3*nbasis*nbasis,
     &                             'vbscf.m','virtual','kscr',
     &                             IGMEM_DEBUG)
         call jaccov(Q(kprma),newbas0,Q(ksao),
     &               Q(kvec),Q(kval),3,crilow,
     &               Q(kscr))
c
c...  harmonic
c
         if (oharm) then
_IF(debug)
            write(iwr,*) '((4)) VIRTUAL - OHARM 2 __ OK'
_ENDIF
            call expharm(Q(kvec),'vectors',ilifq)
            call expharm(flop,'ctrans',flop)
c...     we just make a few zero virtuals (toetoe should sort them out)
            call tdown(Q(kvec),ilifq,Q(kvec),ilifq,nbasis)
c...     end of harmonic dimensions
            otran = .true.
         end if
c
c...  end of *harmonic*
c
         if (iprins.ge.50000) then
            write(iwr,*) ' eigenvectors/values'
            call prvc(Q(kvec),nbasis,nbasis,
     &                Q(kval),'v','l')
         end if
c
         call gmem_free_inf(kscr,'vbscf.m','virtual','kscr')
         call gmem_free_inf(kval,'vbscf.m','virtual','kval')
c
         if (super) then
_IF(debug)
            write(iwr,*) '((5)) VIRTUAL - SUPER TRUE __ OK'
_ENDIF
            nvirt=nbasis - ndoubly
            kvir = kvec + nbasis*(nocc-nsingly)
         else
_IF(debug)
            write(iwr,*) '((6)) VIRTUAL - SUPER FALSE __ OK'
_ENDIF
            nvirt=nbasis-nocc
            kvir = kvec + nbasis*nocc
         end if
c
         call fmove(Q(kvir),virt,nvirt*nbasis)
c
         call gmem_free_inf(kvec,'vbscf.m','virtual','kvec')
         call gmem_free_inf(kprma,'vbscf.m','virtual','kprma')
         call gmem_free_inf(koccp,'vbscf.m','virtual','koccp')
         call gmem_free_inf(ksao,'vbscf.m','virtual','ksao')
c
      else if (aos) then
_IF(debug)
         write(iwr,*) '((13)) VIRTUAL - AOS'
_ENDIF
         if (oharm) call caserr('aos and harmon do not mix')
         call vclr(virt,1,nbasis*nbasis)
         do i=1,nbasis
            virt(i,i)=1.0d0
         end do
         nvirt = nbasis
c
      else if (super_hybrid) then
c
c...  generate the ao's as virtuals as required by super_hybrid
c...  sao and scratch are abused as integers
c
_IF(debug)
           write(iwr,*) '((14)) VIRTUAL - SUPER_HYBRID __ OK'
_ENDIF
           ksh1 = igmem_alloc_inf(nbasis,'vbscf.m','virtual','ksh1',
     &                            IGMEM_DEBUG)
           ksh2 = igmem_alloc_inf(nbasis,'vbscf.m','virtual','ksh2',
     &                            IGMEM_DEBUG)
           if (.not.hybry) 
     &        call caserr('hybrids required for superhybrid')
           reduce = .false.
           kk = 0
c...      note nbasis may not be enough
           call vclr(virt,1,nbasis*nbasis)
           call vclr(Q(ksh1),1,nbasis)
           call vclr(Q(ksh2),1,nbasis)
c
           if (igno_sel.eq.99) then
              if (iprins.ge.10.or.nitscf.le.0) 
     1        write(iwr,'(a)')' *** super hybrid orthog for doubles ***'
c....           generate orthog complement per hybrid
              nvirt = 0
              iatom = 0
9876          iatom = iatom+1
                if (nacat(iatom).gt.0) then
                  if(nvirt+nopa(iatom)-nacat(iatom).gt.maxvirt)
     1            call caserr('marcin-1 - not enough space for virts')
                  kvhyb = igmem_alloc_inf(nopa(iatom)*nopa(iatom),
     &                                    'vbscf.m','virtual',
     &                                    'kvhyb',IGMEM_DEBUG)
                  call virthyb(iatom,occ,virt(1,nvirt+1),
     &                         Q(kvhyb),nopa(iatom),nbasis,
     &                         Q(1))
                  call gmem_free_inf(kvhyb,'vbscf.m','virtual','kvhyb')
                  nvirt = nvirt + nopa(iatom)-nacat(iatom)
                end if
                if (equivir.and.(iacat(1,iatom).eq.eqvmo(1))) then
                  if (nacat(iatom).ne.1) call caserr('equiv nacat ne 1')
                  nn = nopa(iatom)-nacat(iatom)
                  if (nvirt+nn*(neqmo-1).gt.maxvirt)
     1            call caserr('marcin-2 - not enough space for virts')
                  call  primequiv(virt(1,nvirt-nn+1),nbasis,nn,
     1                            'virt','virt')
                  iatom = iatom + neqmo-1
                  nvirt = nvirt + (neqmo-1)*nn
                end if
              if (iatom.lt.natom) go to 9876 
              if (iprins.ge.10.or.nitscf.le.0) 
     1        write(iwr,'(i10,a,i10)') nvirt,' virtuals used - max',
     1                                 maxvirt
c
c...  nopa # ao's ; iopa(j,iatom) index of ao 
c...  nacat # mo' : iacat(j,iatom) index of mo
c...     generate excitation patterns
c...     keep iexw and iex in sync
c...     for now only works for doubles; we handle singles later jvl2008
c
              nn = 0
              kk = nocc
              nequi = 0
              iatom = 0
              iiex = 0
9877          iatom = iatom + 1
                if (nacat(iatom).gt.0) then
                   kstart = kk
                   do imo=1,nacat(iatom)
                      if (.not.ofrzmo(iacat(imo,iatom))) then
                         kk = kstart
                         nexcit = nopa(iatom) - nacat(iatom)
                         nequi = nequi + 1
                         if (eqvmo(1).eq.iacat(1,iatom)) then
c...              equivalent set
                            if (nacat(iatom).ne.1) 
     1                        call caserr('nacat.ne.1 equiv')
                            iequi(nequi) = neqmo
                         else
                            iequi(nequi) = 1
                         end if
c
                         do ieq=1,iequi(nequi)
                            iiex = iiex + 1
                            nex(iiex) = nexcit
                            nexw(iiex) = nexcit
                            iex(1,iiex) = iacat(imo,iatom)
                            iexw(1,iiex) = iacat(imo,iatom)
                            do iao=1,nexcit
                               kk = kk + 1
                               iex(iao+1,iiex) = kk  
                               iexw(iao+1,iiex) = kk
                            end do
                            if (iequi(nequi).gt.1) iatom = iatom + 1
                         end do
                      end if
                   end do
                   if (iequi(nequi).gt.1) iatom = iatom - 1
                end if
              if (iatom.lt.natom) go to 9877 
              if (iprins.ge.10.or.nitscf.le.0)  then
                 write(iwr,9881) nequi,natom
9881             format(/' hand equiv excitations / equic',i6,
     1                  ' atoms',i6)
                 kk = 0
                 do ieq=1,nequi
                    write(iwr,9882) ieq,iequi(ieq)
9882                format(' ****** equiv set',i3,'  of size',i3,
     1                     ' -- excit --')
                    do imo=1,iequi(ieq)
                       kk = kk + 1
                       write(iwr,9883) (iex(k,kk),k=1,nex(kk)+1)
9883                   format(1x,i3,' -- ',15i4,/,(16x,15i4))
                    end do
                 end do
              end if
           else 
c
              do i=1,natom
c...  do not generate for fragments without mo's
                 if (nacat(i).gt.0) then
_IF(debug)
                    write(iwr,*) '((14a)) VIRTUAL - NACAT > 0 __ OK'
_ENDIF
                    if (igno_sel.eq.2) then
c...            sort iopa, to have largest coeff in mo last, accuracy .1
                       do j=nopa(i),nopa(i)-igno_sh(i)+1,-1
                          l = iopa(j,i)
                          dd = 0.0d0
                          do jj=1,nacat(i)
                            jjj = iacat(jj,i)
                            dd = max(dd,abs(occ(l,ncore+jjj)))
                          end do
                          jomax = j
                          do jo=j-1,1,-1
                             ddd = 0.0d0
                             l = iopa(jo,i)
                             do jj=1,nacat(i)
                                jjj = iacat(jj,i)
                                ddd = max(ddd,abs(occ(l,ncore+jjj)))
                             end do
                             if (ddd.gt.dd+0.01d0) jomax = jo
                          end do
                          if (jomax.ne.j) then
c...                  swap
                             l = iopa(jomax,i)
                             iopa(jomax,i) = iopa(j,i)
                             iopa(j,i) = l
                          end if
                       end do
                    end if
c
c...         skip last igno_sh ao's
c
                    do j=1,nopa(i)-igno_sh(i)
                       k = iopa(j,i)
                       if (Q(ksh1+k-1).eq.0.0d0) then
                          Q(ksh1+k-1) = 1.0d0
                          kk = kk + 1
                          virt(k,kk) = 1.0d0
                          Q(ksh2+k-1) = kk
                       end if
                    end do
c...         check if skipped ao's are in mo's
                    write(iwr,601) i,(iopa(j,i),
     1                             j=nopa(i)-igno_sh(i)+1,nopa(i))
601                 format(' ** super hybrid atom ',i4,
     1                     ' - skipped aos ',20i4)
                    do j=nopa(i)-igno_sh(i)+1,nopa(i)
                       k = iopa(j,i)
                       dd = 0.0
                       do jj=1,nacat(i)
                          jjj = iacat(jj,i)
                          dd = max(dd,abs(occ(k,ncore+jjj)))
                       end do
                        if (dd.lt.1.0d-5) 
     1                     call caserr('insufficient ignored ao')
                     end do
                  end if
               end do
               nvirt = kk
c...     generate excitation patterns
c...     keep iexw and iex in syncx
               nn = 0
               do iatom=1,natom
                  do imo=1,nacat(iatom)
                     if (.not.ofrzmo(iacat(imo,iatom))) then
_IF(debug)
                        write(iwr,*) '((14b))VIRTUAL - NOT OFRZMO __ OK'
_ENDIF
                        nn = nn + 1
                        nex(nn) = nopa(iatom) - igno_sh(iatom)
                        nexw(nn) = nopa(iatom) - igno_sh(iatom)
                        iex(1,nn) = iacat(imo,iatom)
                        iexw(1,nn) = iacat(imo,iatom)
                        do iao = 1,nex(nn)
                           kk = Q(ksh2+(iopa(iao,iatom))-1)+0.1d0
                           iex(iao+1,nn) = kk  + nocc
                           iexw(iao+1,nn) = kk  + nocc
                         end do
                      end if
                   end do
                end do
             end if
c
             ncol = nvirt + nocc 
c
             call gmem_free_inf(ksh2,'vbscf.m','virtual','ksh2')
             call gmem_free_inf(ksh1,'vbscf.m','virtual','ksh1')
      end if
c
c...  assign new functions to atoms, if necessary
c
         if (iprins.ge.10000) then
            write(iwr,*) 'orbitals to mix'
            call prvc(virt(1,nbasis-nvirt+1),nvirt,nbasis,virt,'o','l')
         end if
c
         if (oharm) then
_IF(debug)
            write(iwr,*) '((15)) VIRTUAL - HARMONIC __ OK'
_ENDIF
c
c...  move and check zero vectors i.c.o. *harmonic* just to be neat
c
            nn = 0
            do i=nvirt,1,-1
               do j=1,nbasis
                  if (virt(j,i).ne.0.0d0) go to 100
               end do
c...  zero vector found
               nn = nn + 1
               do k=i,nvirt-nn
                  call dcopy(nbasis,virt(1,k+1),1,virt(1,k),1)
               end do
               call vclr(virt(1,nvirt-nn+1),1,nbasis)
100            continue
            end do
            if (nn.ne.newbas1-newbas0)call caserr('vb harmonic problem')
         end if
c
c...   canonicalise doublies
c
       call vbcanon(occ,nbasis,nocc,canonicalise(1),'occ',Q(1))
c
c...   canicalise virtuals
c
       call vbcanon(virt,nbasis,nvirt,canonicalise(2),'virt',
     &              Q(1))
c
_IF(debug)
      write(iwr,*) '((KONIEC)) VIRTUAL'
_ENDIF
c
c...       virt is occ+virt for normal
c
      if (normal) then
         do i=nvirt+nocc,nocc+1,-1
           call fmove(virt(1,i-nocc),virt(1,i),nbasis)
         end do
         call fmove(occ,virt,nbasis*nocc)
      end if
c
      return
      end

      subroutine vbcanon(v,ndim,nnmdim,canonic,vtext,qq)
c
c...  canonicalise the orbital set over the operator icanon
c...  canonic = 0 : nothing
c...  canonic = 1 : v (h-t)
c...  canonic = 2 : fock 
c...  canonic = 10 : localise
c...  canonic = 100 : orthogonalise: doc - Lowdin, voc on doc - Schmidt
c
      implicit REAL (a-h,o-z), integer (i-n)
c
INCLUDE(../m4/common/sizes)
INCLUDE(../m4/common/iofile)
INCLUDE(../m4/common/infoa)
c...  for dynamic memory allocation debugging purpose, contains IGMEM_ vars.
INCLUDE(../m4/common/gmemdata)
INCLUDE(../m4/common/gmempara)
c
INCLUDE(common/turtleparam)
INCLUDE(common/vbcri)
INCLUDE(common/splice)
INCLUDE(common/tractlt)
INCLUDE(common/infato)
INCLUDE(common/vbvirt)
INCLUDE(common/twice)
INCLUDE(common/scftvb)
c
      dimension v(ndim,nnmdim),qq(*)
      integer canonic
      character*(*) vtext
      character*5 toper
c..    trust fortran90 to sort this (* else use mxorbvb *)
      integer ip(nnmdim),ipb(nnmdim)
c
      nmdim = nnmdim
      if (canonic.eq.0) then
         toper = ' '
      else if (canonic.eq.1) then
         toper = 'h'
      else if (canonic.eq.2) then
         toper = 'fock'
      else if (canonic.eq.10) then
         toper = 'local'
      else if (canonic.eq.100) then
         toper = 'ortho'
         if (vtext.ne.'occ') call caserr('wrong set ortho')
      else 
         call caserr(' wrong operator')
      end if
_IF(debug)
        write(iwr,*) '((1)) canon - CANONICALISE 1'
_ENDIF
      ns = ndim * (ndim + 1) / 2
      koper    = igmem_alloc_inf(ns,'vbscf.m','vbcanon','koper',
     &                          IGMEM_DEBUG)
      kvb7_fcao = kscra7vb('kvb7_fcao',ns,'r','n')
      if (toper.eq.' '.or.nmdim.eq.0) then
         call gmem_free_inf(koper,'vbscf.m','vbcanon','koper')
         return
      else if (toper.eq.'h'.or.
     1        (toper.eq.'fock'.and.kvb7_fcao.lt.0)) then
         if (toper.eq.'fock') write(iwr,'(3a)') 
     1   ' *** canon uses h for fock for ',vtext,' this once ***'
         kscr    = igmem_alloc_inf(ns,'vbscf.m','vbcanon','kscr',
     &                             IGMEM_DEBUG)
         call get1e(qq(koper),dummy,'h',qq(kscr))
         call gmem_free_inf(kscr,'vbscf.m','vbcanon','kscr')
      else if (toper.eq.'fock') then
         if (kvb7_fcao.lt.0) call caserr('No fock matrix available')
         kvb7_fcao = kscra7vb('kvb7_fcao',ns,'r','r')
         call rdedx(qq(koper),ns,kvb7_fcao,num8)
      else if (toper.eq.'local'.or.toper.eq.'ortho') then
         call get1e(qq(koper),dummy,'s',qq(koper))
      end if
c
c...  make sure we stay within hybrids
c
      if (hybry) call clmetr(qq(koper),ndim,'any')
c
      kv =  1
      if (vtext.eq.'occ') then
c
c...  for occupied get core and doubly occupied in front
c
         do i=1,nmdim
            ip(i) = 0
         end do
         do i=1,ncore
            ip(i) = i
         end do
         do i=ncore+1,ncore+ndoubly
            ip(i) = idoubly(i-ncore) + ncore
         end do
         khulp = igmem_alloc_inf(ndim,'vbscf.m','vbcanon','khulp',
     1                           IGMEM_DEBUG)
         khul1 = igmem_alloc_inf(nmdim,'vbscf.m','vbcanon','khul1',
     1                           IGMEM_DEBUG)
         khul2 = igmem_alloc_inf(nmdim,'vbscf.m','vbcanon','khul2',
     1                           IGMEM_DEBUG)
         call coperm(ip,qq(khulp),nmdim,iwr)
         call pinver(ip,ipb,nmdim)
         call change(ip,qq(khul1),qq(khul2),
     &               v,1,nmdim,qq(khulp),
     1               ndim,nmdim)
         nmdim = ndoubly
         kv = ncore + 1
      end if
c
      if (toper.eq.'local') then
c...    localise
_IF(debug)
         write(iwr,*) '((10)) CANON - LOCALISE'
_ENDIF
         kvecin = igmem_alloc_inf(ndim*nmdim,'vbscf.m','vbcanon',
     &                            'kvecin',IGMEM_DEBUG)
         kiorot = igmem_alloc_inf(nmdim*nmdim,'vbscf.m','vbcanon',
     &                            'kiorot',IGMEM_DEBUG)
         ktao = igmem_alloc_inf(ndim*(ndim+1)/2,'vbscf.m',
     &                          'vbcanon','ktao',IGMEM__DEBUG)
         ktran = igmem_alloc_inf(nmdim*nmdim,'vbscf.m','vbcanon',
     &                           'ktran',IGMEM_DEBUG)
         krij = igmem_alloc_inf(nat*nmdim*(nmdim+1)/2,'vbscf.m',
     &                          'vbcanon','krij',IGMEM_DEBUG)
         call gen_iorot(qq(kiorot),qq(koper),qq(ktao),v(1,kv),
     1                  nmdim,ndim,hybry)
c
         call fmove(v(1,kv),qq(kvecin),ndim*nmdim)
         call poporb_vb(qq(kvecin),qq(ktran),v(1,kv),
     *                  qq(krij),qq(koper),qq(kiorot), ndim,nmdim,
     *                  (nmdim*nmdim+1)/2,nat,nit_loc)
c
         call gmem_free_inf(krij,'vbscf.m','vbcanon','krij')
         call gmem_free_inf(ktran,'vbscf.m','vbcanon','ktran')
         call gmem_free_inf(ktao,'vbscf.m','vbcanon','ktao')
         call gmem_free_inf(kiorot,'vbscf.m','vbcanon','kioro')
         call gmem_free_inf(kvecin,'vbscf.m','vbcanon','kvecin')
c
      else if (toper.eq.'ortho') then
c      
c...     orthogonalise the doubly occupied orbitals first
c 
        ksmo = igmem_alloc_inf(ndoubly*(ndoubly+1)/2,'vbscf.m',
     &                         'vbcanon','ksmo',IGMEM_DEBUG)
        kcmq = igmem_alloc_inf(ndoubly*ndoubly,'vbscf.m',
     &                         'vbcanon','kcmq',IGMEM_DEBUG)
        klam = igmem_alloc_inf(ndoubly,'vbscf.m',
     &                         'vbcanon','klam',IGMEM_DEBUG)
        kvmt = igmem_alloc_inf(ndim,'vbscf.m',
     &                         'vbcanon','kvmt',IGMEM_DEBUG)
c
        call lowdins(v(1,kv),qq(koper),qq(ksmo),qq(kcmq),
     1               qq(klam),qq(kvmt),ndoubly,ndim,crilow)
c
        call gmem_free_inf(kvmt,'vbscf.m','vbcanon','kvmt')
        call gmem_free_inf(klam,'vbscf.m','vbcanon','klam')
        call gmem_free_inf(kcmq,'vbscf.m','vbcanon','kcmq')
        call gmem_free_inf(ksmo,'vbscf.m','vbcanon','ksmo')
c
c...     orthogonalise the doublies onto the singlies and normalise the singlies
c
        kscr = igmem_alloc_inf(ndim,'vbscf.m',
     &                         'vbcanon','kscr',IGMEM_DEBUG)
c
        call schmids(v(1,ndoubly+kv),v(1,kv),qq(koper),qq(kscr),
     1               nsingly,ndoubly,ndim,crilow,'norm')
c
        call gmem_free_inf(kscr,'vbscf.m','vbcanon','kscr')
        nmdim = ndoubly + nsingly
c
      else
c
c...  standard canonicalisation
c...  get operator in basis of orbitals 
c
         nsm = nmdim*(nmdim+1)/2
         kopm = igmem_alloc_inf(nsm,'vbscf.m','vbcanon','kopm',
     &                          IGMEM_DEBUG)
         kcm = igmem_alloc_inf(ndim,'vbscf.m','vbcanon','kcm',
     &                         IGMEM_DEBUG)
         call maket(qq(kopm),v(1,kv),nmdim,ndim,qq(koper),qq(kcm))
         call gmem_free_inf(kcm,'vbscf.m','vbcanon','kcm')
         if (iprins.ge.50000) then
            write(iwr,*) toper,' operator in basis of ',vtext
            call tripri(qq(kopm),nmdim)
         end if
c
c...  diagonalise the operator
c
         kvec = igmem_alloc_inf(nmdim*nmdim,'vbscf.m','vbcanon','kvec',
     1                          IGMEM_DEBUG)
         kval = igmem_alloc_inf(nmdim,'vbscf.m','vbcanon','kval',
     1                          IGMEM_DEBUG)
         kiky = igmem_alloc_inf(nmdim+1,'vbscf.m','vbcanon','kiky',
     1                          IGMEM_DEBUG)
         kscr = igmem_alloc_inf(nmdim*2,'vbscf.m','vbcanon','kscr',
     1                          IGMEM_DEBUG)
         call filiky(qq(kiky),nmdim+1)
c
         call jacobt(qq(kopm),qq(kiky),nmdim,qq(kvec),
     1                  nmdim,qq(kval),0,2,crilow,qq(kscr))
c
         call gmem_free_inf(kscr,'vbscf.m','vbcanon','kscr')
         call gmem_free_inf(kiky,'vbscf.m','vbcanon','kiky')
         if (iprins.ge.50000) then
            write(iwr,*) toper,' eigenvectors/values'
            call prvc(qq(kvec),nmdim,nmdim,qq(kval),'v','a')
         end if
c
c...  transform the orbitals to the gaussian basis
c
         kvecc = igmem_alloc_inf(ndim*nmdim,'vbscf.m','vbcanon','kvecc',
     1                           IGMEM_DEBUG)
         call fmove(v(1,kv),qq(kvecc),ndim*nmdim)
         call vclr(v(1,kv),1,ndim*nmdim)
         call mxmb(qq(kvecc),1,ndim,qq(kvec),1,nmdim,
     &             v(1,kv),1,ndim,ndim,nmdim,nmdim)
         call gmem_free_inf(kvecc,'vbscf.m','vbcanon','kvecc')
         call gmem_free_inf(kval,'vbscf.m','vbcanon','kval')
         call gmem_free_inf(kvec,'vbscf.m','vbcanon','kvec')
         call gmem_free_inf(kopm,'vbscf.m','vbcanon','kopm')
c
      end if
c
      if (iprins.ge.50000) then
         write(iwr,*) toper,' eigenvectors for ',vtext
         call prvc(v(1,kv),nmdim,ndim,v,'v','a')
      end if
c
      if (vtext.eq.'occ') then
c
c...  Due to canonicalisation over Fock matrix the order of the doubly
c...  occupied orbitals might have changed. Change the order back.
c
         nmdim = nnmdim
         kold = igmem_alloc_inf(ndim*(ncore+nactiv),'vbscf.m',
     1                          'vbcanon','kold',IGMEM_DEBUG)
         ksmo = igmem_alloc_inf((ncore+ndoubly)*(ncore+ndoubly),
     1                          'vbscf.m','vbcanon','ksmo',IGMEM_DEBUG)
         call reorderc(v,ndim,nmdim,qq(kold),qq(ksmo),ip,ncore+ndoubly,
     1                 qq(khul1),qq(khul2),qq(khulp))
         call gmem_free_inf(ksmo,'vbscf.m','vbcanon','ksmo')
         call gmem_free_inf(kold,'vbscf.m','vbcanon','kold')
c
c...  for occupied get doubly occupied in old position
c
         call change(ipb,qq(khul1),qq(khul2),v,1,nmdim,qq(khulp),
     1               ndim,nmdim)
         call gmem_free_inf(khul2,'vbscf.m','vbcanon','khul2')
         call gmem_free_inf(khul1,'vbscf.m','vbcanon','khul1')
         call gmem_free_inf(khulp,'vbscf.m','vbcanon','khulp')
      else if (vtext.eq.'virt')  then
c
c...  got to have phases normalised  for virtuals (virsig)
c
        call virsig(v,ndim,nmdim)
      else
        call caserr(' no occ or virt in vbcanon')
      end if
c
      call gmem_free_inf(koper,'vbscf.m','vbcanon','koper')
c
      return
      end
      subroutine virthyb(iatom,occ,virt,vhyb,nopaa,nbasis)
c
c...   routine to make virtuals for one  hybrid set iatom
c...   no idempotent, do not canonicalise anything
c...  
c...   iatom : the atom we are dealing with
c...   occ   : all occupieds
c...   virt  : the virtuals to be generated
c...   vhyb  : nopaa*nopaa compressed array
c...   nopaa : nopa(iatom) = # ao's in hybrid
c...   nbasis: # ao's in total

c
      implicit none
c
INCLUDE(common/turtleparam)
INCLUDE(../m4/common/sizes)
INCLUDE(common/infato)
INCLUDE(common/vbequiv)
INCLUDE(common/scftvb)
INCLUDE(common/vbcri)
INCLUDE(../m4/common/gmemdata)
INCLUDE(../m4/common/gmempara) 
INCLUDE(../m4/common/vcore) 
INCLUDE(../m4/common/harmon)
INCLUDE(../m4/common/iofile)
c
      integer iatom,j,jj,k,i,nvirt,nbasis,nopaa
      integer itri,igmem_alloc_inf
      integer kshyb,kprma,ksao,kscr,kval
      REAL occ(nbasis,*),virt(nbasis,*),vhyb(nopaa,nopaa)
      REAL dummy
      itri(i,j) = max(i,j)*(max(i,j)-1)/2+min(i,j)
c
c...  nopa # ao's ; iopa(j,iatom) index of ao 
c...  nacat # mo' : iacat(j,iatom) index of mo
c
      if (oharm) call caserr('virthyb not for harmonic yet')
      if (nopaa.ne.nopa(iatom)) call caserr('boem nopaa')
c
      kshyb = igmem_alloc_inf(itri(nopaa,nopaa),'vbscf.m','virthyb',
     &                       'kshyb',IGMEM_DEBUG)
      kprma = igmem_alloc_inf(itri(nopaa,nopaa),'vbscf.m','virhyb',
     &                        'kprma',IGMEM_DEBUG)
c
c...  compress and check mo's
c
      do j=1,nacat(iatom)
         do jj=1,nopaa
            vhyb(jj,j) = occ(iopa(jj,iatom),iacat(j,iatom))
         end do
      end do
c
c...  get metric
c
      ksao = igmem_alloc_inf(itri(nbasis,nbasis),'vbscf.m','virthyb',
     &                       'ksao',IGMEM_DEBUG)
      kscr = igmem_alloc_inf(itri(nbasis,nbasis),'vbscf.m','virthyb',
     &                       'kscr',IGMEM_DEBUG)
c
      call get1e(Q(ksao),dummy,'s',Q(kscr))
      do i=1,nopa(iatom)
         do j=1,i
            Q(kshyb-1+itri(i,j)) = 
     1      Q(ksao-1+itri(iopa(i,iatom),iopa(j,iatom)))
         end do
      end do
c
      call gmem_free_inf(kscr,'vbscf.m','virthyb','kscr')
      call gmem_free_inf(ksao,'vbscf.m','virthyb','ksao')
c      
      call vclr(qq(kprma),1,itri(nopaa,nopaa))
      k=0     
      do i=1,nopaa
         do j=1,i
            k = k + 1 
            do jj=1,nacat(iatom)
               qq(kprma+k-1) = qq(kprma+k-1) + vhyb(i,jj)*vhyb(j,jj)
            end do
         end do
      end do
c
      if (iprins.ge.50000) then
         write(iwr,*) ' projection operator'
         call tripri(qq(kprma),nopaa)
      end if
c
c...  diagonalise the projection operator, virtuals have eigenvalue zero
c
      kval  = igmem_alloc_inf(nopaa,'vbscf.m','virthyb','kval',
     &                        IGMEM_DEBUG)
      kscr    = igmem_alloc_inf(2*itri(nopaa,nopaa)+3*nopaa*nopaa,
     &                          'vbscf.m','virthyb','kscr',IGMEM_DEBUG)
c....    scr = itri(n,n)+itri(n,n)+n*n+2*n*n (about)
      call jaccov(Q(kprma),nopaa,Q(kshyb),
     &            vhyb,Q(kval),3,crilow,
     &            Q(kscr))
c
      if (iprins.ge.50000) then
         write(iwr,*) ' eigenvectors/values'
         call prvc(vhyb,nopaa,nopaa,Q(kval),'v','l')
      end if  
c
c.... expand
c
      nvirt = nopaa - nacat(iatom)
      call vclr(virt,1,nvirt*nbasis)
      do i=1,nvirt
         do j=1,nopaa
            virt(iopa(j,iatom),i)= vhyb(j,i+nacat(iatom))
         end do
      end do
c
         
      call gmem_free_inf(kscr,'vbscf.m','virthyb','kscr')
      call gmem_free_inf(kval,'vbscf.m','virthyb','kval')
      call gmem_free_inf(kprma,'vbscf.m','virthyb','kprma')
      call gmem_free_inf(kshyb,'vbscf.m','virthyb','kshyb')
c
      return
      end
c***********************************************************************
_EXTRACT(poporb,mips4)
**== originally poporb.f
      subroutine poporb_vb(vecin,tran,vecout,rij,sao,iorot,
     *                     l1,n1,n2,natoms,maxit)
c
c     poporb calculates the localized molecular orbitals by the
c     population method of j. pipek and p. g. mezey,
c     j. chem. phys. 90, 4916 (1989),
c
c     modified version for use within vb (jvl 2002)
c     taken from GAMESS (anala); all reads/prints removed
c     the orbitals to be localised are first in vecin/vecout
c
c     sao     : ao s-matrix (l2)
c     vecin   ; input vectors overwritten as well with loc vectors
c     tran    : unitary tranformation matrix to localised orbitals (n1*n1)
c     vecout  : resulting localised vectors
c     rij     : localisation indicators (n2,natoms)
c     l1      : dimension of vectors
c     n1      : # vectors to localise
c     n2      : triangle(n1)
c
c     iir(l1),map(l1),iord(l1),qpix(l1),qpjx(l1) : (formal parameters) moved to /junk/
c     the real arrays are now dynamic arrays (to avoid alignment troube)
c
c     iorot   : allows restriction of rotations (1 : allowed; 0 forbidden)
c               *note* deviates from original norot
c     iorot(n1,n1) :  may be set on basis of 1-electron attraction
c
c
      implicit REAL (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
INCLUDE(common/turtleparam)
INCLUDE(../m4/common/sizes)
INCLUDE(../m4/common/iofile)
INCLUDE(common/vbcri)
      common /junk / limlow(maxat),limsup(maxat),
     1               iir(mxorbvb),map(mxorbvb),iord(mxorbvb)
c
      dimension qpix(mxorbvb),qpjx(mxorbvb)
c
      dimension vecin(l1,n1),tran(n1,n1),sao(*),iorot(n1,n1),
     +          vecout(l1,n1),rij(n2,natoms)
c
c
      parameter (dzero=0.0d+00, done=1.0d+00, two=2.0d+00, four=4.0d+00)
c
      save ifirst
      data ifirst/1/
c
      ind(i,j) = ((max(i,j)*(max(i,j)-1))/2+min(i,j))
c
      tenm10 = cripop
      tenm3 = cripop*1.0d7
      tenm8 = cripop*1.0d2
      cvgloc = cripop*1.0d5
c
      nbasis = l1
      l2 = l1*(l1+1)/2
      l3 = l1*l1
      norb = n1
c
      call dcopy(n1*l1,vecin,1,vecout,1)

c...  determine # basisfunctions per atom (=> common junk)
c
      call aolim
      do 20 k=1,nbasis
      do 20 i=1,natoms
      if(limlow(i).le.k.and.limsup(i).ge.k) map(k)=i
   20 continue
c
      if (norb.eq.0) return
c
      nredo = 0
  110 nredo=nredo+1
c
c          construct initial atomic populations
c
      call vclr(rij,1,n2*natoms)
      ij = 0
      do 280 i = 1,norb
         do 280 j = 1,i
            ij = ij+1
            do 260 k = 1,nbasis
               kk = map(k)
               do 260 l = 1,nbasis
                  ll = map(l)
                  kl = ind(k,l)
                  sum = vecin(k,i)*vecin(l,j)*sao(kl)/two
                  rij(ij,kk) = rij(ij,kk) + sum
                  rij(ij,ll) = rij(ij,ll) + sum
  260       continue
  280 continue
c
c          compute initial localization sum
c
      sumrr = dzero
      do 320 i=1,norb
         if(iorot(i,i).eq.0) go to 320
         ii = ind(i,i)
         do 310 k=1,natoms
            sumrr = sumrr + rij(ii,k)**2
  310    continue
  320 continue
      if (ifirst.eq.1) write(iwr,9010) sumrr
c
c          seed the random function, initialize
c          the localization transformation, etc.
c
      donept = done
      if (nredo.eq.1) dd = rndlmo_vb(donept,vecin,l1,n1)
c
      call vclr(tran,1,norb*norb)
      do 340 i = 1,norb
         tran(i,i) = done
  340 continue
      iter = 0
      shift = atan(donept)
c
c          begin localization cycles
c
  360 continue
      change = dzero
      iter = iter+1
      do 380 i = 1,norb
         iir(i) = i
  380 continue
      nnn = norb
      do 400 i = 1,norb
         dd = rndlmo_vb(change,vecin,l1,n1)
         iii = int(dd*dfloat(nnn)+done)
         iord(i) = iir(iii)
         iir(iii) = iir(nnn)
         nnn = nnn-1
  400 continue
c
c        for each pair of orbitals a two dimensional unitary
c        transformation is performed. the transformation is
c
c           psi'(i) =  cos(t)*psi(i) + sin(t)*psi(j)  and
c           psi'(j) = -sin(t)*psi(i) + cos(t)*psi(j).
c
c        localization requires that t be such as to maximize
c        the sum of the squares of the atomic populations.
c
      do 920 iii = 1,norb
      i = iord(iii)
      if(iorot(i,i).eq.0) go to 920
      ii = ind(i,i)
      jm = 1
      rm = dzero
      tm = dzero
      sm = dzero
      cm = done
      do 580 j = 1,norb
      if(i.eq.j) go to 580
      if(iorot(i,j).eq.0) go to 580
      ij = ind(i,j)
      jj = ind(j,j)
      t = dzero
      tx = dzero
      do 480 kk = 1,natoms
         t = t + four*rij(ij,kk)**2 - rij(ii,kk)**2 - rij(jj,kk)**2
     *         + two*rij(ii,kk)*rij(jj,kk)
         tx = tx + rij(ij,kk)*(rij(jj,kk) - rij(ii,kk))
  480 continue
      if ((dabs(t) .le. tenm10) .and. (dabs(tx) .le. tenm10)) go to 580
      tx = four*tx
      t = datan2(tx,t)/four
      sign = done
      if (t .gt. dzero) sign = -done
      t = t+sign*shift
      itim = 0
  500 itim = itim+1
      s = dsin(t)
      cc = dcos(t)
      rin = dzero
      do 520 kk = 1,natoms
         qpi = cc*cc*rij(ii,kk)+s*s*rij(jj,kk)+two*cc*s*rij(ij,kk)
         qpj = cc*cc*rij(jj,kk)+s*s*rij(ii,kk)-two*cc*s*rij(ij,kk)
         rin = rin+qpi*qpi+qpj*qpj-rij(ii,kk)**2-rij(jj,kk)**2
  520 continue
      ttest = dabs(t)-shift
      if ((dabs(t) .le. tenm8) .or. (dabs(ttest) .le. tenm8)) go to 560
      if (rin .ge. -tenm8) go to 560
      if (itim .le. 1) go to 540
         write (iwr,9020) i,j
         write (iwr,9030) t,s,cc,rin
      return
c
  540 sign = done
      if (t .gt. dzero) sign = -done
      t = t+shift*sign
      go to 500
c
  560 if (rin .le. rm) go to 580
      rm = rin
      tm = t
      sm = s
      cm = cc
      jm = j
  580 continue
c
      rin = rm
      t = tm
      s = sm
      cc = cm
      j = jm
      ij = ind(i,j)
      jj = ind(j,j)
      if(iorot(i,j).eq.0) go to 920
c
c        accumulate the 2x2 rotation
c
      change = change+t*t
      call drot(norb,tran(1,i),1,tran(1,j),1,cc,s)
c
c        update the atomic populations
c
      do 880 kk = 1,natoms
         qpi = cc*cc*rij(ii,kk)+s*s*rij(jj,kk)+two*cc*s*rij(ij,kk)
         qpj = cc*cc*rij(jj,kk)+s*s*rij(ii,kk)-two*cc*s*rij(ij,kk)
         qpij = (cc*cc-s*s)*rij(ij,kk)+cc*s*(rij(jj,kk)-rij(ii,kk))
         do 720 k = 1,norb
            if (i.eq.k.or.j.eq.k) goto 720
            ik = ind(i,k)
            jk = ind(j,k)
            qpix(k) = cc*rij(ik,kk)+s*rij(jk,kk)
            qpjx(k) = cc*rij(jk,kk)-s*rij(ik,kk)
            rij(ik,kk) = qpix(k)
            rij(jk,kk) = qpjx(k)
  720    continue
         rin = rin+qpi+qpj-rij(ii,kk)-rij(jj,kk)
         rij(ii,kk) = qpi
         rij(jj,kk) = qpj
         rij(ij,kk) = qpij
  880 continue
  920 continue
c
c          test for convergence of localization procedure
c
      change = sqrt(two*change/(norb*(norb-1)))
      if(ifirst.eq.1.and.mod(iter,20).eq.0) write(iwr,9050) iter,change
      If(iter.lt.maxit  .and.  change.gt.tenm3*cvgloc) go to 360
      if(change.le.cvgloc) go to 1000
c        if(nredo.le.2) write(iwr,9060)
c        if(nredo.le.2) go to 110
            if (ifirst.eq.1) write(iwr,9070)
            go to 1010
c
c          finished with localization cycles
c
 1000 continue
      if (ifirst.eq.1) write (iwr,9080) iter,change
c
c        transform to final orbitals, copy virtual space
c
1010  continue
      do 930 i = 1 , norb
         call vclr(rij(1,1),1,nbasis)
         do 940 j = 1 , norb
         sc = tran(j,i)
         call daxpy(nbasis,sc,vecout(1,j),1,rij(1,1),1)
 940     continue
         call dcopy(nbasis,rij(1,1),1,vecin(1,i),1)
 930  continue
c
c ----- now load up lmos into vecin array
c
      do 270 i = 1 , norb
         call dcopy(nbasis,vecin(1,i),1,vecout(1,i),1)
 270  continue
c ----
c     if(oprint(54)) then
c        call prsq(tran,norb,norb,n1)
c     end if
c
c          construct final atomic populations
c
      call vclr(rij,1,n2*natoms)
      ij = 0
      do 1280 i = 1,norb
         do 1280 j = 1,i
            ij = ij+1
            do 1260 k = 1,nbasis
               kk = map(k)
               do 1260 l = 1,nbasis
                  ll = map(l)
                  kl = ind(k,l)
                  sum = vecout(k,i)*vecout(l,j)*sao(kl)/two
                  rij(ij,kk) = rij(ij,kk) + sum
                  rij(ij,ll) = rij(ij,ll) + sum
 1260       continue
 1280 continue
c
c               compute final localization sum (skip frozen mo-s)
c
      sumrr = dzero
      do 1380 i=1,norb
         if(iorot(i,i).eq.0) go to 1380
         ii = ind(i,i)
         do 1370 k = 1,natoms
            sumrr = sumrr + rij(ii,k)**2
 1370    continue
 1380 continue
      if (ifirst.eq.1) write(iwr,9140) sumrr
      ifirst = 0
c
      return
c
 9010 format(/10x,'the initial localization sum is',f14.6)
 9020 format(1x,'no rotation increases atomic populations',
     *           ' --- localization aborted'/
     *           10x,'i=',i3,5x,'j=',i3)
 9030 format(5x,8htheta = ,g20.10/5x,12hsin(theta)= ,f10.7,
     +     15h   cos(theta)= ,f10.7/5x,29htotal change to this point = ,
     +     g20.10)
 9050 format(10x,'iteration',i5,'   orbital change=',g20.10)
 9060 format (//10x,'*** localization has been unsucessful ***'
     +        //10x,'program will restart with new random number'
     +         /10x,'and rotation sequence for orbitals')
 9070 format( 10x,'+++++++++++++++++++++++++++++++++++++++++++++++'/
     +        10x,'+           localization stopped              +'/
     +        10x,'+++++++++++++++++++++++++++++++++++++++++++++++')
 9080 format(' localization converged in',i5,' it. with change ',E11.2)
 9140 format(10x,'the final localization sum is',f16.6/)
      end
**==rndlmo.f

c***********************************************************************
      function rndlmo_vb(dxx,d,l1,na)
c  
c     mod to get rid of infoa dependence (na)
c     + changed some funny suff with n,m (jvl2002 )
c
      implicit REAL (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
c  
      dimension d(l1,na),u(1)
c  
      save u  
c  
      parameter (dzero=0.0d+00, done=1.0d+00)
c  
      pi = acos(-done)
      if (dxx .ne. dzero) then
         n = abs(na-l1)+1
         m = mod(n+5,na) + 1 
         dxy = d(n,m)*atan(done)
         u(1) = (pi+dxy)**5
         dxy = dfloat(int(u(1)))
         u(1) = u(1)-dxy
         rndlmo_vb = u(1)
         return  
      endif   
c  
      u(1) = (pi+u(1))**5
      dxy = dfloat(int(u(1)))
      u(1) = u(1)-dxy
      rndlmo_vb = u(1)
      return  
      end     

c***********************************************************************
      subroutine gen_iorot(iorot,hao,tao,virt,nvirt,nbasis,hybry)
c
c...  generate a iorot for poporb_vb that does not mix e.g. sigma-pi
c...  we might have to include the occupieds to generste better classes
c
      implicit REAL (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
c
      logical hybry
      dimension hao(*),tao(*),virt(nbasis,nvirt),iorot(nvirt,nvirt)
INCLUDE(common/turtleparam)
INCLUDE(../m4/common/sizes)
INCLUDE(common/vbcri)
      common/junk/ vv(mxorbvb)
c
      call get1e(tao,dummy,'v',hao)
      if (hybry) call clmetr(tao,nbasis,'any')
      call maket(hao,virt,nvirt,nbasis,tao,vv)
c
      ij = 0
      do i=1,nvirt
        do j=1,i
          ij = ij + 1
          if (dabs(hao(ij)).gt.critor) then
             iorot(i,j) = 1
             iorot(j,i) = 1
          else
             iorot(i,j) = 0
             iorot(j,i) = 0
          end if
        end do
      end do
c
      return
      end

c***********************************************************************
_ENDEXTRACT
      subroutine fockvb(q)
c
      implicit REAL (a-h,o-z),integer(i-n)
c
c...  This subroutine builds a Fock matrix on ao basis (including
c...  virtuals.
c
      dimension q(*)
INCLUDE(common/scftvb)
INCLUDE(../m4/common/sizes)
INCLUDE(../m4/common/iofile)
INCLUDE(../m4/common/files)
INCLUDE(common/ffile)
INCLUDE(common/turtleparam)
INCLUDE(common/c8_16vb)
      dimension iky(melec)
INCLUDE(common/splice)
INCLUDE(common/tractlt)
INCLUDE(../m4/common/gmemdata)
INCLUDE(../m4/common/gmempara) 
      common/blkin/gg(340),mij(340),mword,mdumm
_IF(linux)
      external fget
_ENDIF
c             t   
c      korneel   = igmem_alloc_inf(nbasis*(nbasis+1)/2,'vbscf.m','vbscf',
c     &                       'ks',IGMEM_DEBUG)
c...  Gather information to partition core memory
      call wr15vb(nelec,nalfa,nconf,ndet,nstruc,maxdet,maxstr,
     &            nwpack,ncoeff,imax,ndum11,ndum12,'read')
      kdens1   = igmem_alloc_Inf(lenbas,'vbscf.m','fockvb',
     &                           'kdens1',IGMEM_DEBUG)
      kfao = igmem_alloc_Inf(lenbas,'vbscf.m','fockvb',
     &                           'kfao',IGMEM_DEBUG)
      kh = igmem_alloc_Inf(lenbas,'vbscf.m','fockvb',
     &                           'kh',IGMEM_DEBUG)
      kscr = igmem_alloc_Inf(lenbas,'vbscf.m','fockvb',
     &                           'kscr',IGMEM_DEBUG)
c...  read one electrondensity matrix on ao basis
      kvb7_dn = kscra7vb('kvb7_dn',lenbas,'r','n')
      if (kvb7_dn.lt.0) 
     1   call caserr('1-el densmatrix needed in fockvb but not found')
      kvb7_dn = kscra7vb('kvb7_dn',lenbas,'r','r')
      call rdedx(q(kdens1),lenbas,kvb7_dn,num8)
c...  Create Fock matrix on ao basis
      call vclr(q(kfao),1,lenbas)
c...  get H-core
      call get1e(q(kh),dummy,'h',q(kscr))
      if (iprint.gt.1500) then
        write(iwr,*) 'one electron h-matrix'
        call tripri(q(kfao),nbasis)
      end if
c...  Adjust the 1-electron density matrix
      k=0
      do i=1,nbasis
         do j=1,i
            if (i.eq.j) then
              q(kdens1+k)=q(kdens1+k)/4
            else
              q(kdens1+k)=q(kdens1+k)/2
            endif
            k=k+1
         enddo
      enddo
c...  Start building Fock matrix on ao basis by running over
c...  two electron integrals
        do 20113 i=1,n2file
           lbl=n2blk(i)
           if (lbl.eq.n2last(i)) go to 20113
           iunit=n2tape(i)
           call search(n2blk(i),iunit)
68          call fget(gg,k,iunit)
           if (k.ne.0) then
              call sgmatvb(q(kfao),q(kdens1))
              lbl=lbl+1
              if (lbl.ne.n2last(i)) go to 68
           end if
20113   continue
_IF(parallel)
         call pg_dgop(24,q(kfao),lenbas,'+')
_ENDIF
         call vadd(q(kfao),1,q(kh),1,q(kfao),1,lenbas)
c
      if (iprint.gt.1500) then
         write(iwr,*) 'fock matrix on ao basis'
         call tripri(q(kfao),nbasis)
      end if
c...  Store fock matrix on ao basis on ed7
      kvb7_fcao = kscra7vb('kvb7_fcao',lenbas,'r','w')
      call wrt3(q(kfao),lenbas,kvb7_fcao,num8)
c
      call gmem_free_set(kdens1,kscr)
c
      return
      end

      subroutine overlon(nelec,nalfa,igroup,ngroup,idet,q,lword)
c
      implicit REAL (a-h,o-z), integer (i-n)
c
c...  This subroutine calculates the overlap of the "old" wavefunction
c...  with the new one
c
      dimension q(*),igroup(5,*),idet(*)
c
INCLUDE(../m4/common/iofile)
INCLUDE(common/tractlt)
INCLUDE(common/splice)
INCLUDE(../m4/common/sizes)
INCLUDE(common/turtleparam)
INCLUDE(common/brill)
INCLUDE(common/scftvb)
      kbcoef = 1
      ndet = igroup(1,1)
_IF(atmol)
      rewind 33
c...  read "new" coefficients
      read(33) ndum,(q(kbcoef),ij=1,ndet)
_ELSE
      kvb7_bcn = kscra7vb('kvb7_bcn',igroup(1,1),'r','n')
      if (kvb7_bcn.lt.0) then
         call caserr('Psi0 coefs needed in ovelon but not found')
      else
         kvb7_bcn = kscra7vb('kvb7_bcn',igroup(1,1),'r','r')
         call rdedx(q(kbcoef),igroup(1,1),kvb7_bcn,num8)
      endif
_ENDIF
      if (nitscf.gt.1) then
c
c...     determine the overlap of the current wave function with the
c...     previous one
c
         kbold = kbcoef + ndet
c
c...     read old coefficients
c
_IF(atmol)
c...     read "old" coefficients
         read(33) (q(ij),ij=kbold,ndet)
_ELSE
         kvb7_bco = kscra7vb('kvb7_bco',igroup(1,1),'r','n')
         if (kvb7_bco.lt.0) then
            call caserr('Old Psi0 coefs needed in ovelon but not found')
         else
            kvb7_bco = kscra7vb('kvb7_bco',igroup(1,1),'r','r')
            call rdedx(q(kbold),igroup(1,1),kvb7_bco,num8)
         endif
_ENDIF
c
         kvold = kbold + ndet
c
c.....   read old and new vectors
c
_IF(atmol)
         read(33) nbasis,nv,(q(ij),ij=kvold,kvold+nbasis*nv-1)
         if (nv.ne.ncore+nscf) call caserr('tape33 goof')
_ELSE
         kvb7_vo = kscra7vb('kvb7_vo',nbasis*(ncore+nscf),'r','n')
         if (kvb7_vo.lt.0) then
            call caserr('Old vectors needed in overlon but not found')
         else
            kvb7_vo = kscra7vb('kvb7_vo',nbasis*(ncore+nscf),'r','r')
            call rdedx(q(kvold),nbasis*(ncore+nscf),kvb7_vo,num8)
         endif
_ENDIF
         kvnew = kvold + nbasis*(ncore+nscf)
_IF(atmol)
         read(33) nbasis,nv,(q(ij),ij=kvnew,kvnew+nbasis*nv-1)
_ELSE
         kvb7_vn = kscra7vb('kvb7_vn',nbasis*(ncore+nscf),'r','n')
         if (kvb7_vn.lt.0) then
            call caserr('New vectors needed in overlon but not found')
         else
            kvb7_vn = kscra7vb('kvb7_vn',nbasis*(ncore+nscf),'r','r')
            call rdedx(q(kvnew),nbasis*(ncore+nscf),kvb7_vn,num8)
         endif
_ENDIF
c
         ksmo = kvnew + nbasis*(ncore+nscf)
         ksao = ksmo + nscf*nscf
c        ksao = ksmo + nbasis*nbasis
         kscr = ksao + nbasis*nbasis
         ktot = kscr + 2*max(nbasis,nscf)*nbasis
c        ktot = kscr + 2*nbasis*nbasis
c
c...     this is somewhat pessimistic
c
         if (ktot.gt.lword) call corfait(ktot,lword,'in brigen/1')
c
         call get1e(q(ksao),dummy,'s',q(ksao))
c
c...     skip core orbitals
c
         kvold = kvold + ncore*nbasis
         kvnew = kvnew + ncore*nbasis
c
c...     old with itself
c
         call oldnew(q(kbold),q(kbold),ndet,idet,
     &            nelec,nalfa,q(kvold),q(kvold),q(ksao),q(ksmo),q(kscr),
     &            nbasis,nscf,swavaa)
c
c....    new with itself
c
         call oldnew(q(kbcoef),q(kbcoef),ndet,idet,
     &            nelec,nalfa,q(kvnew),q(kvnew),q(ksao),q(ksmo),q(kscr),
     &            nbasis,nscf,swavbb)
c
c....    old with new
c
         call oldnew(q(kbold),q(kbcoef),ndet,idet,
     &            nelec,nalfa,q(kvold),q(kvnew),q(ksao),q(ksmo),q(kscr),
     &            nbasis,nscf,swavab)
c
c....    normalise (might not be necessary)
c
         swave = swavab / (dsqrt(swavaa) * dsqrt(swavbb))
      end if
c
      return
      end

c***********************************************************************
      subroutine virsig(v,n,nmdim)
c
c...  This subroutine looks for the largest coefficient of 
c...  a virtual orbital (in 2 decimals) and makes sure that this
c...  coefficient becomes a positive number. 
c
      implicit none
      integer n,i,j,maxcol,nmdim
      REAL v(n,nmdim)
c
      do i=1,nmdim
         maxcol=0
         do j=1,n
            if (int(100*dabs(v(j,i))).gt.abs(maxcol)) then
               maxcol=int(100*(v(j,i)))
            endif
         enddo
         if (maxcol.lt.0) then
            do j=1,n
               v(j,i)=-1.0d0*v(j,i)
            enddo
         endif
      enddo
      return
      end

c***********************************************************************
      subroutine makepsi0(v,ci)
c
c...  This subroutine freezes the orbitals that are occupied twice
c...  always. Then it transforms the integrals from the raw (ao) basis
c...  to the mo basis. After that it generates the H and S matrices on
c...  structure basis. Then it calculates the ci-vector or psi(0). Finally
c...  it calculates the ao density matrix (to be used as convergence 
c...  criterium).
c...  Parameters :
c...    v     - vectors (mo's)
c...    lword - total number of words available to makepsi0
c...    ci    - vb vector
c...    qq    - dynamicly allocated vectors
      implicit none
c
c... for dynamic memory allocation debugging purpose, contains IGMEM_ vars.
INCLUDE(../m4/common/gmemdata)
INCLUDE(../m4/common/gmempara)
INCLUDE(../m4/common/vcore)
c
INCLUDE(../m4/common/atmol3)
c
INCLUDE(../m4/common/sizes)
INCLUDE(../m4/common/timeperiods)
INCLUDE(../m4/common/iofile)
INCLUDE(../m4/common/runlab)
INCLUDE(common/turtleparam)
INCLUDE(common/c8_16vb)
INCLUDE(common/tractlt)
INCLUDE(common/splice)
INCLUDE(common/hsmattype)
INCLUDE(common/brill)
INCLUDE(common/scftvb)
INCLUDE(common/basisvb)
INCLUDE(common/ffile)
INCLUDE(common/vbtimess)
      common /txema/ struc,ivec
      logical struc
      integer ivec,ijj,ijk,isbass,nwpi
c
      integer igmem_alloc_inf,igmem_alloc_all_inf,igmem_max_memory
      integer kscra7vb
      integer nelec,nalfa,nconf,ndets,nstruc,maxdet,maxstr,
     &        nwpack,ncoeff,imax,manual,ndoub,oldncore,oldnsa,
     &        oldncurt,ngroup,northo,nn,ktemp,kvec
      integer kh,ks,kg,kscr,kpacd,kndet,kidps,kcoef,
     &        kigro,kiort,kidet,kjdet,kjjdet,kipos,kdao,kdmo,kno,kdete,
     &        ksao,khao,kvnat,kdoub,ksing,kindex,kidoub
      integer oldncol,i,ndettot,lensec,oldprint,itemp
      integer kmaxmem,kmemscr,kmemleft,nipw
      integer kvb7_transvb
      character*8 oldrunt
      REAL oldcore
      REAL  v(*),a1,a0,cpulft,ci(*)
c...  ci - ci vector passed from parent subroutine 

      if (scfconv) then
         if (zruntp.ne.'scf') then
            if (p0core) then
               print *,' ** frozen core density needs fixing **' 
               p0core = .false.
            end if
         end if
c...     upon convergence do not include doublies in core to make print better
         p0core = .false.
      end if
c
c...  print info
c
      if (iprint.ge.10) then
         write(iwr,*) 'generating psi(0)...'
      endif
c
c...  read symbolic information from crestr
c
      call wr15vb(nelec,nalfa,nconf,ndets,nstruc,maxdet,maxstr,
     &            nwpack,ncoeff,imax,manual,ndoub,'read')
c
c...  keep certain parameters 
c
      oldcore = core
      oldrunt = zruntp
      oldncol = ncol
      oldncore = ncore
      oldnsa = nsa
      oldncurt = ncurt
      zruntp = 'scf'
      if (p0core) then
         hsmatt = 'vbfunction'
         ncurt = -1
         ncore = ncore + ndoub
      endif
      ncol = oldncore + nscf
      nsa = ncol - ncore
      struc=.true.
      lenact = nsa*(nsa+1)/2
*
      if (unitary) then
         call orthovb(v,nbasis,nsa,ncore)
      end if
*
c
c...  if not done reserve space on ed7 for vectors (we have oldnsa and nsa)
c...  used in transformvb -- super hybrid can be unpredictable
c
c...     reserve nbasis*nbasis +nsa*nbasis, because this might be necessary later !
      kvb7_transvb = kscra7vb('kvb7_transvb',
     1                        (nbasis+max0(nsa,oldnsa))*nbasis,'r','w')
c
      kh     = igmem_alloc_inf(lenact,'vbscf.m','makepsi0','kh',
     &                         IGMEM_DEBUG)
      kvec   = igmem_alloc_inf(ncol*nbasis,'vbscf.m','makepsi0',
     &                         'kvec',IGMEM_DEBUG)
      kidoub = igmem_alloc_inf(ndoub,'vbscf.m','makepsi0','kidoub',
     &                         IGMEM_DEBUG)
_IF(debug)
      write(iwr,*) '[(1)] makepsi0' 
_ENDIF
c
      call fmove(v,Q(kvec),ncol*nbasis)
c
      ks     = igmem_alloc_inf(lenact,'vbscf.m','makepsi0','ks',
     &                         IGMEM_DEBUG)

_IF(debug)
      write(iwr,*) '[(5)] makepsi0' 
_ENDIF
c
c...  check if we need to reorder the orbitals to put doubles between core and active
c
      if (p0core.and.ndoub.gt.0) then
_IF(debug)
      write(iwr,*) '[(55)] makepsi0' 
_ENDIF
         call swapvec(v,Q(kvec),Q(kidoub),
     &                oldncore,ndoub,ncol,nbasis)
      endif
c
c...  transform the integrals (1st transformation only variably occupied orbitals)
c
      a1 = cpulft(1)
      call start_time_period(TP_VB_TRAN)

      call transformvb(Q(ks),Q(kh),
     &                 Q(kvec))

      kg     = igmem_alloc_inf(n2int+1,'vbscf.m',
     &                         'makepsi0','kg',IGMEM_DEBUG)

      call getin2(Q(kg))
      call clredx


_IF(debug)
      write(iwr,*) '[(6)] makepsi0' 
_ENDIF
c
c...  print integrals over orbitals if requested
c   
      if (iprint.ge.50) then
         write(iwr,*) ' 1-electron integrals over orbitals'
         call tripri(Q(kh),nsa)
         write(iwr,*) ' overlap matrix between orbitals'
         call tripri(Q(ks),nsa)
         if (iprint.ge.100000) then
            write(iwr,*) ' 2-electron integrals'
            call tripri(Q(kg+1),lenact)
         end if
      end if
      a1 = cpulft(1) - a1
      t4indx0 = t4indx0 + a1
      n4indx0 = n4indx0 + 1
      call end_time_period(TP_VB_TRAN)
c
c...  build h/s-matrix for psi(0)
c
      a0 = cpulft(1)
      call start_time_period(TP_VB_ME)

      call hsmat(Q(ks),Q(kh),Q(kg),nsa,Q(1))
      a1 = cpulft(1)
      tmatre0 = tmatre0 + a1 - a0
      nmatre0 = nmatre0 + 1
      call end_time_period(TP_VB_ME)
_IF(debug)
      write(iwr,*) '[(7)] makepsi0' 
_ENDIF
      kmaxmem = igmem_max_memory()
      kmemscr = kmaxmem/10
      kmemleft = kmaxmem - kmemscr
      kscr = igmem_alloc_inf(kmemscr,'vbscf.m','makepsi0','kscr',
     &                       IGMEM_DEBUG)
c
c...  diagonalise h-matrix
c
      a1 = cpulft(1)
      call start_time_period(TP_VB_DIAG)

c...  DAVIDT SUBROUTINE STILL USING OLD MEMORY PARTITIONING
      call davidt(Q(kscr),kmemscr,core,evb,'psi0')
c
c...  move ci-vector to right position
c
      do i=1,nstruc
         ci(i)=Q(kscr+i-1)
      end do

      a1 = cpulft(1) - a1
      tdavid = tdavid + a1
      call end_time_period(TP_VB_DIAG)
_IF(debug)
      write(iwr,*) '[(8)] makepsi0' 
_ENDIF
      call gmem_free_inf(kscr,'vbscf.m','makepsi0','kscr')
c
      kidet = igmem_alloc_inf(nwpi(nelec*ndets),'vbscf.m',
     &                        'makepsi0','kidet',IGMEM_DEBUG)
      kjdet = igmem_alloc_inf(nwpi(nelec*ndets),'vbscf.m',
     &                        'makepsi0','kjdet',IGMEM_DEBUG)
      kjjdet = igmem_alloc_inf(nwpi(nelec),'vbscf.m',
     &                        'makepsi0','kjjdet',IGMEM_DEBUG)
      kpacd = igmem_alloc_inf(nwpack,'vbscf.m','makepsi0','kpacd',
     &                         IGMEM_DEBUG)
      kndet = igmem_alloc_inf(nstruc,'vbscf.m','makepsi0','kndet',
     &                        IGMEM_DEBUG)
      kidps = igmem_alloc_inf(ndets*nstruc,'vbscf.m','makepsi0','kidps',
     &                        IGMEM_DEBUG)
      kcoef = igmem_alloc_inf(ncoeff,'vbscf.m','makepsi0','kcoef',
     &                        IGMEM_DEBUG)
      kigro = igmem_alloc_inf(5*nstruc,'vbscf.m','makepsi0','kigro',
     &                        IGMEM_DEBUG)
      kiort = igmem_alloc_inf(ncore+nsa,'vbscf.m','makepsi0','kiort',
     &                        IGMEM_DEBUG)
      kmaxmem = igmem_max_memory()
      kmemscr = kmaxmem/10
      kscr  = igmem_alloc_inf(kmemscr,'vbscf.m','makepsi0','kscr',
     &                            IGMEM_DEBUG)
c
c...  find out about orthogonality classes
c
      call sym1(Q(ks),nsa,Q(kiort),northo,'noprint')
_IF(debug)
      write(iwr,*) '[(9)] makepsi0' 
_ENDIF
c   remco
c...  read information from datatape
c...  INFORM SUBROUTINE STILL USING OLD MEMORY PARTITIONING
c
      call inform(nelec,nalfa,ndets,nstruc,Q(kpacd),
     &            Q(kndet),Q(kidps),Q(kcoef),
     &            ncoeff,Q(kigro),ngroup,Q(kiort),
     &            northo,Q(kscr),kmemscr)
      call gmem_free_inf(kscr,'vbscf.m','makepsi0','kscr')

_IF(debug)
      write(iwr,*) '[(10)] makepsi0' 
_ENDIF
c
c...  Perform PCM tweak option (if specified)
c
      if (ntscf.eq.1.and.ipcmt.gt.0) then
_IF(debug)
      write(iwr,*) '[(10a)] makepsi0' 
_ENDIF
        if (nitscf.eq.0) then
          itemp = maxscf_s
          maxscf_s = maxscf
          maxscf = itemp
        end if
c...  FIXED - originaly q() instead ci()
        if (ipcmt.lt.1000) then
          do i=1,nstruc
            if (i.eq.ipcmt) then
              ci(i)=sqrt(0.9)
            else
              ci(i)=sqrt(0.1/(nstruc-1))
            endif
          enddo
        endif
        write(iwr,1000) ipcmt
1000  format('WARNING: VB-vector is tweaked. Mainly structure',i3,
     &   ' is present this VBSCF run!')
      end if
      if (ntscf.eq.2.and.ipcmt.gt.0.and.nitscf.eq.0) then
_IF(debug)
      write(iwr,*) '[(10b)] makepsi0' 
_ENDIF
        itemp = maxscf
        maxscf = maxscf_s
        maxscf_s = itemp
      endif
_IF(debug)
      write(iwr,*) '[(11)] makepsi0' 
_ENDIF
c
c...  generate psi0 on determinant basis (after this qq(kscr) contains
c...  the coefficients in the psi0 on det. basis. ndettot is the total
c...  number of determinants).

c...  PSI0DET SUBROUTINE STILL USING OLD MEMORY PARTITIONING
c
      kdete = igmem_alloc_inf(2*ndets+nelec,'vbscf.m','makepsi0',
     &                        'kdete',IGMEM_DEBUG)

      call psi0det(ci,nstruc,Q(kigro),ngroup,Q(kcoef),
     &             ncoeff,ndets,Q(kndet),Q(kidps),
     &             Q(kidet),Q(kjdet),Q(kpacd),
     &             nelec,nalfa,ndettot,Q(kdete),
     &             kmemleft,'save') 

_IF(debug)
      write(iwr,*) '[(12)] makepsi0' 
_ENDIF
c
c...  move determinant coefficient immediately after the ci-vector
c 
      kdmo  = igmem_alloc_inf(nsa*(nsa+1)/2,'vbscf.m',
     &                        'makepsi0','kdmo',IGMEM_DEBUG)
      kdao  = igmem_alloc_inf(max(nbasis,ncore+nsa)*nbasis,
     &                        'vbscf.m','makepsi0','kdao',IGMEM_DEBUG)
      kipos = igmem_alloc_inf(nelec*nelec,'vbscf.m','makepsi0','kipos',
     &                        IGMEM_DEBUG)
      kvnat = igmem_alloc_inf(max(ncol,nbasis)*nbasis,
     &                        'vbscf.m','makepsi0','kvnat',IGMEM_DEBUG)
c
c...  unpack information of mo's per determinant
c
      call unpack(Q(kpacd),n8_16,Q(kidet),ndettot*nelec)
_IF(debug)
      write(iwr,*) '[(13)] makepsi0' 
_ENDIF
c
c...  If converged make sure to run the natorbt in the correct mode:
c...  In scf mode only the 1-electrondensity matrix is calculated, in
c...  geometry optimisation mode also the 2-electrondensity matrix and
c...  the lagrangian are needed.
c
      if (scfconv) then
         zruntp = oldrunt
      endif
c
c...  first two-electron integral is (00/00) !!!
c...  qq(kg) => qq(kg+1)
c
      call natorbt(Q(kdete),ndettot,Q(kidet),
     &             Q(kjdet),Q(kjjdet),nelec,nalfa,
     &             Q(kdao),Q(kdmo),Q(kipos),
     &             Q(kvnat),Q(kiort),northo,nbasis,nsa,
     &             ncore,Q(kh),Q(ks),Q(kg+1))

_IF(debug)
      write(iwr,*) '[(14)] makepsi0' 
_ENDIF

c
c...  move original parameters back
c
      struc=.false.
      zruntp = oldrunt
      hsmatt = 'full      '
      ncore = oldncore
      ncol = oldncol
      nsa = oldnsa
      lenact = nsa*(nsa+1)/2
      ncurt = oldncurt
      core_i = core
      core = oldcore

c...  free whole bunch of memory - one after another for debugging purpose
      call gmem_free_inf(kvnat,'vbscf.m','makepsi0','kvnat')
      call gmem_free_inf(kipos,'vbscf.m','makepsi0','kipos')
      call gmem_free_inf(kdao,'vbscf.m','makepsi0','kdao')
      call gmem_free_inf(kdmo,'vbscf.m','makepsi0','kdmo')
      call gmem_free_inf(kdete,'vbscf.m','makepsi0','kdete')
      call gmem_free_inf(kiort,'vbscf.m','makepsi0','kiort')
      call gmem_free_inf(kigro,'vbscf.m','makepsi0','kigro')
      call gmem_free_inf(kcoef,'vbscf.m','makepsi0','kcoef')
      call gmem_free_inf(kidps,'vbscf.m','makepsi0','kidps')
      call gmem_free_inf(kndet,'vbscf.m','makepsi0','kndet')
      call gmem_free_inf(kpacd,'vbscf.m','makepsi0','kpacd')
      call gmem_free_inf(kjjdet,'vbscf.m','makepsi0','kjjdet')
      call gmem_free_inf(kjdet,'vbscf.m','makepsi0','kjdet')
      call gmem_free_inf(kidet,'vbscf.m','makepsi0','kidet')
      call gmem_free_inf(kg,'vbscf.m','makepsi0','kg')
      call gmem_free_inf(ks,'vbscf.m','makepsi0','ks')
      call gmem_free_inf(kidoub,'vbscf.m','makepsi0','kidoub')
      call gmem_free_inf(kvec,'vbscf.m','makepsi0','kvec')
      call gmem_free_inf(kh,'vbscf.m','makepsi0','kh')

_IF(debug)
      write(iwr,*) '[(KONIEC)] makepsi0' 
_ENDIF

      return
      end

c***********************************************************************
      subroutine swapvec(vin,vout,idoub,ncore,ndoub,nscf,nbasis)
c
      implicit none
c
INCLUDE(../m4/common/iofile)
      integer off1,off2,i,j,k,ndoub,idoub(ndoub),ncore,nbasis,nscf
      integer kscra7vb,k7doub
      REAL vin(nbasis*nbasis),vout(nbasis*nbasis)
      if (ndoub.eq.0) call caserr('swapvec unsuitable for ndoub.eq.0')
      k7doub = kscra7vb('k7doub',ndoub,'i','r')
      call readi(idoub,ndoub,k7doub,num8)
      call fmove(vin,vout,ncore*nbasis)
      j=1
      k=1
      do i=1,nscf - ncore
         off1 = (ncore+i-1)*nbasis+1
         if (idoub(j).eq.i) then
            off2 = (ncore+j-1)*nbasis+1
            j=j+1
         else
            off2 = (ncore+ndoub+k-1)*nbasis+1
            k=k+1
         endif
         call fmove(vin(off1),vout(off2),nbasis)
      enddo
      end

c***********************************************************************
      subroutine reorderc(vnew,nbasis,norb,vold,smo,ip,ndoub,h1,h2,hp)
c
c...
c...  When the doubly occupied orbitals are canonicalised over the Fock operator
c...  the order in which they appear might change. DIIS gets confused when the orbital
c...  order is changed or when the phase is changed. This subroutine reorders the
c...  the doubly occupied orbitals and makes the phase equal to previous iterations.
c...
c
      implicit none
INCLUDE(common/turtleparam)
INCLUDE(../m4/common/sizes)
      integer nbasis,norb,i,j,k,l,ip(norb)
      REAL vnew(nbasis,norb),vold(nbasis,norb)
      REAL smo(*),h1(*),h2(*),hp(*)
      REAL smax,signm,ddot
      integer nelec,ndum,ndoub
      integer kscra7vb,kvb7_vo
INCLUDE(../m4/common/iofile)
c
c...  read symbolic information from crestr
c
      kvb7_vo = kscra7vb('kvb7_vo',nbasis*norb,'r','n')
      if (kvb7_vo.lt.0) then
c        first iteration: no old vectors yet
         return
      else
         kvb7_vo = kscra7vb('kvb7_vo',nbasis*norb,'r','r')
         call rdedx(vold,nbasis*norb,kvb7_vo,num8)
      endif
c
c...  get the old doubles together
c
      call change(ip,h1,h2,vold,1,nbasis,hp,nbasis,norb)
      do i=1,ndoub
        do j=1,ndoub
           k=(i-1)*ndoub+j
           smo(k)=ddot(nbasis,vold(1,i),1,vnew(1,j),1)
        enddo
      enddo
      do i=1,ndoub
        smax=-1000000
        do j=1,ndoub
          k=(i-1)*ndoub+j
          if (abs(smo(k)).gt.smax) then
            signm=smo(k)/abs(smo(k))
            smax=abs(smo(k))
            l=j
          end if
        enddo
        do k=1,nbasis
           vold(k,i)=signm*vnew(k,l)
        enddo
      enddo
      call fmove(vold,vnew,ndoub*nbasis)
      return
      end
c
      subroutine primequiv(v,nbasis,norb,set,label)
c
      implicit none
c
      integer nbasis,norb,kmo,iorb 
      character*(*) label,set
      REAL v(nbasis,*)
c
c...   reapply equivalences ; **not nice**
c
INCLUDE(../m4/common/iofile)
INCLUDE(common/vbequiv)
c
      integer iao,imo,mmo,mao
      REAL vmax,vv
      if (.not.equiv) call caserr('equiv confusion')
c
      if (set.eq.'occ') then
         vmax = 0.0d0
c
         do imo=2,neqmo
            do iao=1,neqao
               vv = abs(v(eqvao(iao,imo),eqvmo(imo))
     1                 -v(eqvao(iao,1),eqvmo(1)))
               if (vmax.lt.vv) then
                  vmax = vv
                  mmo = eqvmo(imo)
                  mao = eqvao(iao,imo)
               end if
               v(eqvao(iao,imo),eqvmo(imo)) = v(eqvao(iao,1),eqvmo(1))
            end do
         end do
c
c        write(iwr,1) vmax,mmo,mao,label
         if (vmax.gt.0.0) write(iwr,1) vmax,mmo,mao,label
1        format(/' ** equiv deviation of ',1pe12.5,' at mo',i3,' ao',i3,
     1          5x,a)
c
      else if (set.eq.'virt') then
c...     generate all euivalent virtuals, original set has norb
c...     equivalent oritals have to be consequetive
         kmo = norb 
         call vclr(v(1,kmo+1),1,norb*(neqmo-1)*nbasis)
         do imo=2,neqmo
            if (eqvmo(imo).ne.eqvmo(imo-1)+1) 
     1      call caserr('equiv not in order in primequiv')
            do iorb=1,norb
               kmo = kmo + 1
               do iao=1,neqao
                  v(eqvao(iao,imo),kmo) = v(eqvao(iao,1),iorb)
               end do
            end do
         end do
      else
         call caserr('unrecognised set in primequiv')
      end if
c
      return
      end
*
************************************************************************
      subroutine orthovb(v,nbasis,nmo,ncore)
************************************************************************
*
*     orthogonalise all for unitary
*
      implicit REAL (a-h,o-z), integer (i-n)
*
INCLUDE(../m4/common/gmemdata)
INCLUDE(../m4/common/gmempara)
INCLUDE(../m4/common/vcore)
INCLUDE(common/vbcri)
*
      dimension v(nbasis,*)
*
      nocc = ncore + nmo
      ksao = igmem_alloc_inf(nbasis*(nbasis+1)/2,'vbscf.m','orthovb',
     +                          'ksao',IGMEM_DEBUG)
      kscr = igmem_alloc_inf(nbasis*(nbasis+1)/2,'vbscf.m','orthovb',
     +                          'kscr',IGMEM_DEBUG)
      call get1e(Q(ksao),dummy,'s',Q(kscr))
      call gmem_free_inf(kscr,'vbscf.m','orthovb','kscr')
*
      kscr    = igmem_alloc_inf(nbasis*2,'vbscf.m','orthovb','kscr',
     +                                IGMEM_DEBUG)
*
      call normvc(v,Q(ksao),Q(kscr),nbasis,nocc,cridep)
      call gmem_free_inf(kscr,'vbscf.m','orthovb','kscr')
      call gmem_free_inf(ksao,'vbscf.m','orthovb','ksao')
*
      return
      end
*
************************************************************************
