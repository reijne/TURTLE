
c  $Author: jvl $
c  $Date: 2014-09-15 17:31:14 +0200 (Mon, 15 Sep 2014) $
c  $Locker:  $
c  $Revision: 6296 $
c  $Source: /c/qcg/cvs/psh/GAMESS-UK/vb/vbin.m,v $
c  $State: Exp $
c  $Log: vbin.m,v $
c  Revision 1.89  2007/10/18 16:45:34  jvl
c  Problem with clearing nodes on Huygens fixed (inside subroutine vberr).
c  Some cosmetic changes.
c  Hopefully everything works :)
c  /marcin
c
c  Revision 1.86  2007/09/25 10:35:51  jvl
c  Fixed pg_sndrcv (call to sendrecv which calls MPI_Sendrecv) because of
c  integer*8/integer*4 type mismatch (for big endians), which lead to
c  checksum error.
c  /marcin
c
c  Revision 1.85  2007/08/19 19:54:47  mrdj
c  Use long common block 'sector' first (gfortran tends to truncate to size of
c  first occurrance)
c
c  Revision 1.84  2007/08/19 13:52:12  jvl
c  added level as synonym for shift, as I can't remember
c
c  Revision 1.83  2007/07/06 18:45:44  jvl
c  Minor bug fix in CRESTR with STRUC NORM.
c  Minor update of VBMO directive.
c  Addition of new stripblanks subroutine to strip all the front and trailing
c  blank spaces from string, while keeping all the within.
c  /marcin
c
c  Revision 1.82  2007/07/04 12:37:56  jvl
c  Introducing new VB switch - VBMO, first step to separate output file for molden.
c  Now, upon request, separate file with VB orbitals (with proper format) can be created,
c  called out.vbmo.<title_from_VBMO_switch>. Check vbin.m comment/GAMESS manual for more.
c  Major bugfixes of dimensions in makepsi0 subroutine:
c   - serial + CORE XX
c   - parallel + bigger basis set than 6-31g + CORE switch
c  Minor changes in igmem_alloc_inf names.
c  /marcin
c
c  Revision 1.81  2007/03/20 14:49:31  jvl
c  Pretty major overhaul of the VB code
c  say 50% is now dynamic using igmem_alloc and, oh wonder,
c  all the examples are still checking out OK, including a new resonating
c  super hybrid. (will follow)
c
c  Revision 1.80  2007/03/15 14:03:09  jvl
c  update super hybrid oprion in vb (ignore ao per atom)
c  added damping
c  clean up
c
c  Revision 1.79  2007/02/20 15:21:16  jvl
c  The array atomao (hybrid per ao) made the hybrid code a lot simpler.
c  Unfortunately now we need the possibility to let an ao belong to two hybrids.
c  Thus packing was in order and a special function (ochatao) was made to handle atomao
c  Further some restrictions were removed; now resomating BLW may be used (:-))
c
c  Revision 1.78  2006/12/13 21:56:08  rwa
c  First check in of the experimental VB response property code.
c
c  Revision 1.77  2006/04/27 19:26:58  jvl
c  Upon finding that the orthogonalisations in VB are inconsistent
c  a cleanup was started, which is no way fnished.
c  natorb is doing all sort of non related stuff and makepsi0 calls natorb
c  which is stupid, also frozen core gradients, which are rather trivial
c  for VB should be properly added, etc. but we have to do some work again
c
c  Revision 1.76  2006/04/18 16:35:34  jvl
c  - added basis option to allow Vb to be run in anothr basis.
c  - cleaned up various bits of VB, like core partitioning
c  - adapted get1e
c
c  Revision 1.75  2006/01/13 17:56:48  jmht
c  Merging in changes from the release branch
c
c  Revision 1.74  2005/09/23 14:34:56  jmht
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
c  Revision 1.73.2.2  2005/07/19 06:52:20  wab
c
c  min0/max0 changes (see m4 checkin)
c
c  Revision 1.73.2.1  2005/06/24 05:42:35  wab
c  Second attempt to do this after the bug in guess that caused apparent
c  failures to be diagnosed yesterday has been fixed
c
c  Revision 1.73  2005/05/11 22:37:32  jvl
c  made paralleol vb - I8 completer + few checks and better errormessages
c  sendrecv introduced to avoid race conditions (added to ga-stuff really)
c
c  Revision 1.72  2005/04/22 11:35:04  jvl
c  stats_vb i8 change
c
c  Revision 1.71  2005/04/22 11:07:53  jvl
c  All vb criteria are now in common vbcri and may be set
c  in the input. Also their defaults are now consistent
c
c  Revision 1.70  2005/03/31 14:59:58  jvl
c  - corrected pert option, which was screwed earlier
c  - include file hsinfo
c
c  Revision 1.69  2005/03/07 14:08:07  jvl
c  - made sure that hybrids and core work => change small to any
c    (in old small calculations cores would be OK for hybrids (within small)
c    now for bigger calculations this is not the case anymore
c  - tightened critor and orthogonalisation criteria. 10**-10 is not good enough
c    anymore => 10**-14 or 15
c
c  Revision 1.68  2005/02/05 18:04:25  wab
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
c  Revision 1.67  2005/01/13 14:46:36  rwa
c  Check on vbscf stopcriterion and davidson stop criterion, so that the
c  latter is always smaller than the former.
c
c  Revision 1.66  2005/01/05 10:22:52  jvl
c  Bugfix concerning erasing scratch information between consequetive VBSCF
c  runs (for PCM). Further extension of the PCM model with a tweak input
c  option (JJE).
c
c  Revision 1.65  2004/08/19 15:52:33  jvl
c  Bugfix in virtual fock space. Bugfix in fock canonicalise option. Better
c  ebi print out. Parameter change (JJE)
c
c  Revision 1.64  2004/05/19 22:57:44  jvl
c  Brought harmonic to VB; If hamonic is specified the cartesion
c  components are projected out of the occupied and virtual orbitals
c  yielding a harmonic vbscf and vbci (the latter may of may not be what you want)
c
c  Revision 1.63  2004/05/11 10:27:12  jvl
c  Changed maximum iterations default to 50. Changed default convergence
c  criterium for Brillouin matrix element to 3.0d-5. Changed output of
c  energy to zeroes if maximum number of iterations is reached. Changed
c  maximum number of aos per hybrid from 100 to 150.
c
c  Revision 1.62  2004/03/23 23:54:14  jvl
c  permutations prove too difficult ; now ok i hope
c
c  Revision 1.61  2004/03/23 15:07:55  jvl
c  enhancved conf manual and added default sections
c
c  Revision 1.60  2004/01/29 16:21:38  jvl
c  bugfix in convergence criterion specification (JJE).
c
c  Revision 1.59  2004/01/29 12:58:35  jvl
c  Extended TURTLE with real brillouin criterion, which is now default. Fixed
c  some bugs regarding perturbation option. Changed print-out of iterations:
c  brm value added. (JJE)
c
c  Revision 1.58  2003/11/09 14:53:40  jvl
c  the layout of ed7 in vb was partly redone every time in an geometry optimisation
c  This cause errors and a possible oveflow of ed7 ....
c
c  Revision 1.57  2003/10/22 15:21:43  jvl
c  Fixed a bug concerning "rotten metric" message in case of use of d-type
c  orbitals. Thanks Remco. JJE 22/10/2003
c
c  Revision 1.56  2003/08/28 11:51:04  jvl
c  Tiny bugfix in Makefile.in regarding ranlib. Addition of compilation time
c  in print out.
c
c  Revision 1.55  2003/08/26 13:35:57  jvl
c  Weights are printed together with their branching or rumer diagram.
c  (JJE 26/8/2003)
c
c  Revision 1.54  2003/06/10 16:43:53  jvl
c  Addition of option p0core (JJE 10/6/2003).
c
c  Revision 1.53  2003/04/14 05:50:33  jvl
c  Changed output of vbscf, so that it resembles the original output. Natural
c  orbitals arfe no longer printed by default, only the occupation numbers
c  are printed instead (jje 14/4/2003)
c
c  Revision 1.52  2003/04/07 11:47:43  hvd
c  The WATCOM compiler does not allow goto's into another branch of an
c  if-statement. I have moved the offending call to vberr out of the if-
c  statement until after the goto 1 statement.
c
c  Revision 1.51  2003/02/28 10:58:32  jvl
c  Fixed bug in subroutine makepsi0. Lowered vbscf.m deck from -O3 to -O2
c  for sgi systems (JJE).
c
c  Revision 1.50  2003/02/18 17:17:16  jvl
c  A lot of changes, in general:
c  - In vbci the call to bsmat for natural orbitals is no longer necessary
c  - The suboptions pert and fock for optimise in scf are implemented
c  - The subroutine vbscf is completely restructured and is more logical
c  - Addition of the makepsi0, which uses a small extra 2-electron,
c    transformation
c    JJE - 18 feb 2003
c
c  Revision 1.49  2003/01/19 16:23:52  jvl
c  attempt to correct a comment (in the previous cvs checking) that was interpreted
c  by m4 .....; hope a fixed it ......
c
c  Revision 1.48  2003/01/19 15:52:34  jvl
c  corrected ifdef parallel
c  vb now end properly for hcpu
c
c  Revision 1.47  2003/01/19 12:33:06  jvl
c  - check on # active orbitals (an infinite loop could result)
c  - hcpu option to allow primitive restarting of matrix element calculation in vbci
c
c  Revision 1.46  2002/11/22 14:06:00  jvl
c  reada => rchar / format => format_servec (clash with mopac common)
c  cleaning vb-input ...
c
c  Revision 1.45  2002/11/18 17:35:27  jvl
c  Modifications to add servec, in VB
c  putqvb is now just interface to putq, rdistr requires maximum dimension
c  schmids has 'norm' keyword to denote, if a normalisation is required
c  and few other routines made more straightforward.
c
c  Revision 1.44  2002/11/08 16:01:03  jvl
c  The use of omrdci to clear arrays in unpack has unfortunate consequences
c  for vb; The runtype is set to ci and the mrdci is called when vb is meant
c  so new variable oclunp introduced
c
c  Revision 1.43  2002/11/05 22:43:13  jvl
c  In turtle sometimes 1-electron integrals are calculated before ispchk is called.
c  Therefore this call included .....
c
c  Revision 1.42  2002/10/15 15:15:49  jvl
c  blksiz had changed, so now properly included
c  alse enlarged VB dimensions
c
c  Revision 1.41  2002/09/11 12:07:35  jvl
c  1-electron ints are needed in input for TURTLE, so they are calculated from vbstart
c  if they are not there yet. Also added smallness criterium to atom_hybrid
c
c  Revision 1.40  2002/09/05 14:37:13  jvl
c  1. First modification steps on input for old "pert" option, now opti.
c  2. First modification steps on input for core optimalisation.
c
c  Revision 1.39  2002/09/04 17:01:43  jvl
c  Fixed bug concerning logical guesao. guesao was set by subroutine vbin, but
c  not stored in common block. Afterwards it was used by vbfirst and hybmak.
c  guesao is now stored in common block infato. (JJE 4-9-2002)
c
c  Revision 1.38  2002/08/15 13:20:21  jvl
c  *-1.0d0 is a syntax error according to IBM; fixed
c  outstanding rndlmo_vb fix
c
c  Revision 1.37  2002/05/29 11:24:17  jvl
c  corrected hybrid inputting (nscf => nsa)
c
c  Revision 1.36  2002/05/28 15:07:48  jvl
c  General cleanup
c  got rid of all fortran-files (now part of ed7 is used)
c  handle 3 more common's through include
c  make vector sections in vb more logical
c  added print parameter to getqvb
c  ....
c
c  Revision 1.35  2002/05/20 21:36:42  jvl
c  added a (really long overdue) way to specify hybrids by atom
c
c  Revision 1.34  2002/05/16 15:20:46  jvl
c  added extra directive to combine
c
c  Revision 1.33  2002/05/03 15:53:19  jvl
c  Added screen option to allow cleaning of vectors (in vectors combine)
c
c  Revision 1.32  2002/05/03 12:26:01  jvl
c  Cleaned up getqvb/putqvb nouw adapted vectors are ok as inpu
c  added extra vectors combine option
c  made #its in virtual localisation input parameter
c  poporb does not compile well with high optimisation on sgi
c
c  Revision 1.31  2002/04/28 22:30:28  jvl
c  Added virtual options for vb-orbital optimisation
c  in the process an adapted version of poporb is made
c
c  Revision 1.30  2002/04/25 17:21:33  jvl
c  added checks on negative norm on vector in davidson (= overcomplete sywstem)
c  changed a few (print)flags
c
c  Revision 1.29  2002/02/10 21:07:00  jvl
c  cleaned up printing and made comparison ci/scf possible
c
c  Revision 1.28  2002/02/07 16:14:37  jvl
c  forgot the gop when getting h/s matrices ; fixed
c
c  Revision 1.27  2002/02/07 12:33:32  jvl
c  fixed bugs + added print of orthog h-matrix
c
c  Revision 1.26  2002/01/16 17:02:20  jvl
c  Two clashes of VB with mopac resolved:
c  reada was in bothe : changed in vb ro rchar
c  common/bonds/ in vb was a routine in mopac : changed to common/vbbonds/
c
c  Revision 1.25  2001/12/17 17:28:50  jvl
c  few extra bugfixes found when porting to linuxppc
c
c  Revision 1.24  2001/10/26 11:42:10  jvl
c  bugfix of core directive;
c
c  Revision 1.23  2001/10/24 21:45:39  jvl
c  corrected few problems in geometry optimisation with VB and added initial stuff
c  to limit output later..
c
c  Revision 1.22  2001/10/19 16:30:23  jvl
c  clean hybrid (default) now tries to make hybrid exact
c  root broadcasts orbitals + excitation info for parallel jobs
c  hybrids are checked for consistency
c  infato now contains fragment/ao (atomao) numbers making lots of routines easier
c  some routines had name-changes (clao=> clvecvb)
c  eliminated confusing excitation prints
c
c  Revision 1.21  2001/10/15 15:32:40  jvl
c  major parallel overhaul; makes asynchronous switchable
c  no luck yet however ...........
c
c  Revision 1.20  2001/10/05 07:42:56  jvl
c  stupid error in vb fixed
c
c  Revision 1.19  2001/10/04 12:09:48  jvl
c  added more vector combination options
c
c  Revision 1.18  2001/09/19 16:02:44  jvl
c  some corrections resulting from DEC compile
c
c  Revision 1.17  2001/09/07 14:28:10  jvl
c  common stats in vb clashed with mpi; => stats_vb
c  so about 10 minutes after installing gamess on sara, the first bugfix
c
c  Revision 1.16  2001/08/20 21:14:31  jvl
c  addedd a vbit to vectors combine options (begins to look like servec....)
c
c  Revision 1.15  2001/07/04 21:21:12  jvl
c  various cosmetic improvements; extended super hybrid option
c  p.s. mchines seems to think all files have changed ...
c
c  Revision 1.14  2001/07/03 15:44:01  jvl
c  6 => iwr
c
c  Revision 1.13  2001/07/03 14:13:16  jvl
c  fixed a few bugs in parallel code and allowed non-usage of peigss in vb
c  (i.e. peigss specifies if peigas is to be used in vb)
c
c  Revision 1.12  2001/07/02 09:16:04  jvl
c  added ignore/exclude to super hybrid
c  added frozen core contributions properly
c  bit of printing fixes
c
c  Revision 1.11  2001/06/30 13:15:48  jvl
c  made vb input bit clearer + bugfixes
c
c  Revision 1.10  2001/06/28 09:06:45  jvl
c  Jacobi davidson checkin voor VB (generalised eigenvalue problem)
c  Xiny 2001
c
c  Revision 1.9  2001/06/27 15:28:47  jvl
c  Changed vector printing
c  fixed various bugs (related to getin2)
c  added frozen en frozen hybrids options (experimental)
c
c  Revision 1.8  2001/06/22 12:28:40  jvl
c  Within the hybrid option it is possible to leave out one or more aos
c  from all hybrid definitions. One could want to do this on purpose. If
c  one leaves out one or more aos an extra hybrid definition, carrying
c  the name *rest* is created. This hybrid has no mo, it only contains
c  all unused aos. (JJE 22jun01)
c
c  Revision 1.7  2001/06/21 16:38:20  jvl
c  added super hybrid option (experimental and perhaps useless) and a bit of cleaning up
c
c  Revision 1.6  2001/06/18 12:52:38  jvl
c  fixed bug + improved vectors print
c
c  Revision 1.5  2001/06/15 15:16:15  jvl
c  added extra in veccomb (vb) + cleanup
c
c  Revision 1.4  2001/06/12 12:21:17  jvl
c  - Removed several bugs
c  - Improved parallel 4-index transformation including changes in getin2
c  - Added timing for VB routines
c
c  Revision 1.3  2001/06/04 21:59:20  jvl
c  added vectors colmbine option + cleaned vb print
c
c  Revision 1.2  2001/04/14 19:59:01  jvl
c  added ? as comment (as it is in gamess)
c
c  Revision 1.1  2001/02/08 14:26:07  jvl
c  Adapted for parallelisation in GAMESS-UK with MPI and GA
c  vbindens.m was splitted into vbin.m and vbdens.m to separate
c  input and density matrix part
c
c  Revision 1.2  2000/03/31 14:18:26  jvl
c  Changed vnorm to vnrm
c  Changed get1e to use getmat
c  Changed scopy to icopy
c  Improved on information about titles
c  Fokke Dijkstra
c
c  Revision 1.1  2000/02/16 12:20:44  jvl
c  adding vb files to gamess repository
c
c Revision 1.16  1998/02/25  13:05:37  fokke
c added get1e for gamess6, new keyword gamess5 uses gamess5 dumpfile
c
c Revision 1.15  1998/02/16  15:14:27  fokke
c removed ifdef MPI
c
c Revision 1.14  1998/02/16  14:49:00  fokke
c included default select 250 when using evec
c
c Revision 1.13  1998/02/16  14:22:44  fokke
c added counter process to parallel implementation
c added output of energies and interactions of single determinants
c
c Revision 1.12  1997/05/30  11:10:21  fokke
c Changed icopy to icopy and removed variables icopy
c
c Revision 1.11  1997/05/23  08:34:25  fokke
c Removed some problems for CPP
c
c Revision 1.10  1997/05/22  12:53:09  joop
c changed imove to icopy
c
c Revision 1.9  1997/05/22  11:31:48  joop
c Added include macro.cpp
c
c Revision 1.8  1997/05/21  11:36:03  joop
c some fixes (21-5-97) FD
c
c Revision 1.7  1997/01/02  16:21:22  joop
c MPI changes for SGI + 1 bugfix
c
c Revision 1.6  1996/11/07  17:10:01  joop
c SiGr => SGI + bugfixje
c
c Revision 1.5  1996/10/29  15:24:01  joop
c rcs log + parallel + labeled matrix elements
c
c
      subroutine vbstart(qq,iwhich)
c
      implicit REAL (a-h,o-z), integer (i-n)
INCLUDE(../m4/common/iofile)
INCLUDE(../m4/common/timeperiods)
INCLUDE(../m4/common/scra7)
c... for dynamic memory allocation debugging purpose, contains IGMEM_ vars.
INCLUDE(../m4/common/gmemdata)
INCLUDE(../m4/common/gmempara)
c
INCLUDE(../m4/common/dump3)
INCLUDE(common/turtleparam)
INCLUDE(common/vbproper)
c
      dimension qq(*)
c...   print LOGO only once
      data idone/0/
      save idone
      character*10 cdate,cname
      character*5 ctime
c     
c...  calculate 1-electron ints if they are not present
c...  may be needed for determining hybrids
c
      call sectst(ionsec,if1)
      if (if1.eq.0) then
         ibasis = ispchk()
         call standv(0,qq)
      end if
c
      call start_time_period(TP_VB)
c
      if (idone.eq.0) then
         write(iwr,11)
_IF(parallel)
11    format(//,40x,50('_'),//
     &,38x,'     t                      t       l              ',/
     &,38x,'     t                      t       l              ',/
     &,38x,'    tttt    u  u   rrrrr   tttt     l         eee  ',/
     &,38x,'     t      u  u    r   r   t       l        eeeee ',/
     &,38x,'     t   t  u  u    r       t   t   l   l    e     ',/
     &,38x,'      ttt    uu     r        ttt     lll      eee  '
     &,///
     &,9x, '  ----              ----       ',
     &     '     ----       ))))) --         ----   )))))      ',/
     &,9x, '       ----                   ',
     &     '----    _ _ _   { o o }  ----   _ _ _   { o o }    ',/
     &,9x, '                     ----     ',
     &     '      _(_)_(_)_  ( ^ ) ---    _(_)_(_)_  ( ^ )     ',/
     &,9x, ' ----         ----             ',
     &     '   _(_)_(_)_(_)_/ /        _(_)_(_)_(_)_/ /        ',/
     &,9x, '         ----              ----',
     &     ' _(_)_(_)_(_)_(_)/       _(_)_(_)_(_)_(_)/         ',/
     &,9x, '                               ',
     &     '(_)_(_)_(_)_(_)_(_)     (_)_(_)_(_)_(_)_(_)        ',/
     &,9x, ' ----     ----      ----       ',
     &     '  //  //   //  //         //  //   //  //           ',
     &          //40x,50('_'))
_ELSE
11    format(//,40x,50('_'),//
     &,38x,'     t                      t       l              ',/
     &,38x,'     t                      t       l              ',/
     &,38x,'    tttt    u  u   rrrrr   tttt     l         eee  ',/
     &,38x,'     t      u  u    r   r   t       l        eeeee ',/
     &,38x,'     t   t  u  u    r       t   t   l   l    e     ',/
     &,38x,'      ttt    uu     r        ttt     lll      eee  '
     &,///
     &,    '  ----              ----       ',
     &     '     ----         ----           ----   )))))      ',/
     &,    '       ----                    ',
     &     '----       ----         ----   _ _ _   { o o }     ',/
     &,    '                     ----      ',
     &     '                   ----      _(_)_(_)_  ( ^ )      ',/
     &,    ' ----         ----             ',
     &     '      ----                 _(_)_(_)_(_)_/ /        ',/
     &,    '         ----              ----',
     &     '                         _(_)_(_)_(_)_(_)/         ',/
     &,    '                               ',
     &     '  ----                  (_)_(_)_(_)_(_)_(_)        ',/
     &,    ' ----     ----      ----       ',
     &     '           ----           //  //   //  //           ',
     &          //40x,50('_'))
_ENDIF
         call vbversion(cdate,ctime,cname)
         write (iwr,12) cdate,ctime,cname
12    format(//,' This version is compiled on ',a10,' ',a5,' by ',a10)
         iflop = kscra7vb('init',ibl7la,'o','i')
      end if
c

c
      idone = idone + 1
c

      if (iwhich.eq.0) then
         kmaxmem = igmem_max_memory() 
         kmemscr = kmaxmem/4
         kmemscr = max(kmemscr,25000)
         kmemscr = min(kmemscr,2000000)
         kmemleft = kmaxmem - kmemscr 
c
         ibase = igmem_alloc_inf(kmemleft,'vbin.m','vbstart','ibase',
     &                           IGMEM_DEBUG) ! MARCIN
         if (kmemleft.lt.50000) then
            write(iwr,101) kmemleft
         else
            write(iwr,102) kmemleft
         end if
101      format (/' WARNING: low on memory, only ',I5,
     &           ' words available for CRESTR',/)
102      format (/' memory available for CRESTR: ',I12,' words'/)
         call crestr(qq(ibase),kmemleft,qq(ibase),qq)
         call gmem_free_inf(ibase,'vbin.m','vbstart','ibase')
      else if (iwhich.eq.1) then
         iflop = kscra7vb('reset',ibl7la,'o','j')
         kmaxmem = igmem_max_memory() 
         kmemscr = kmaxmem/2
         kmemscr = max(kmemscr,25000)
         kmemscr = min(kmemscr,2000000)
         kmemleft = kmaxmem - kmemscr 
c
         ibase = igmem_alloc_inf(kmemleft,'vbin.m','vbstart','ibase',
     &                           IGMEM_DEBUG) ! MARCIN
         if (kmemleft.lt.50000) then
            write(iwr,111) kmemleft
         else
            write(iwr,112) kmemleft
         end if
111      format (/' WARNING: low on memory, only ',I5,
     &           ' words available for WBIN',/)
112      format (/' memory available for VBIN: ',I12,' words'/)
         call vbin(qq(ibase),kmemleft,qq)
         call gmem_free_inf(ibase,'vbin.m','vbstart','ibase')
      else
         iflop = kscra7vb('reset',ibl7la,'o','j')
         kmaxmem = igmem_max_memory() 
c
         if (kmaxmem.lt.50000) then
            write(iwr,121) kmaxmem
         else
            write(iwr,122) kmaxmem
         end if
121      format (/' WARNING: low on memory, only ',I5,
     &           ' words available for VB',/)
122      format (/' memory available for VB: ',I12,' words'/)
         call vb(qq)
      end if
      ibl7la = kscra7vb('end',0,'o','k')
c
      if (iwhich.gt.1.and.oprop) then
         iflop = kscra7vb('reset',ibl7la,'o','j')
         call vbpropdriver(qq)
         ibl7la = kscra7vb('end',0,'o','k')
      endif
      call end_time_period(TP_VB)
c 
      return
      end
      subroutine vb(qq)
c
c...  non-orthogonal valence bond program
c...  ci and scf
c...  version 1.0 1988 (august)
c...  version 2.0 1995
c
c...  J. Verbeek, J. H. Langenberg, C.P. Byrman, F. Dijkstra,
c...  R. W. A. Havenith, J. J. Engelberts, M. L. Zielinski,
c...  Z. Rashid and J. H. van Lenthe, 1988-2013
c...  Theoretical Chemistry Group
c...  Utrecht University
c...  Padualaan 14
c...  3584 CH Utrecht
c...  The Netherlands
c...  Fax: +31-30-2537504
c
      implicit REAL (a-h,o-z), integer (i-n)
      dimension qq(*)
c
c...  adjust ncorq to fit the size of your problem or machine
c
INCLUDE(../m4/common/sizes)
INCLUDE(common/turtleparam)
INCLUDE(common/vbqc)
INCLUDE(common/basisvb)
c     common // q(ncorq)
c
c... for dynamic memory allocation debugging purpose, contains IGMEM_ vars.
INCLUDE(../m4/common/gmemdata)
INCLUDE(../m4/common/gmempara)
c
      parameter (n35=35)
      common /stats_vb/ noccurr(n35)
INCLUDE(common/tractlt)
c
      common /bypass/ index4,ihmat,idavid
c
INCLUDE(common/hsinfo)
c
      common /title/ titel,chtime(4)
      character*80 titel
      character*8 chtime,chyear
c
INCLUDE(common/brill)
c
      common/crestn/ struct
      character*44 struct,dummy
c
INCLUDE(common/vbtimess)
INCLUDE(common/vbpert)
INCLUDE(../m4/common/iofile)
INCLUDE(common/first_vb)
INCLUDE(common/vbproper)
      common/hcpu_vb/iprhs,irowhs
      character*10 charwall
_IF(parallel)
c
c...  this must be first statement of program
c
c     call poe_init()
c     call barrier
_ENDIF
c
c*** TO BE LOOKED AT
c
_IF(parallel)
c     call poe_next
_ENDIF
c
c...  two scratch-files are used : 33 for vectors (length :
c...                               +/- 2*nbasis**2
c...  unit 33 is now replaced by ed7 (num8) with addressing in kscra7vb
c...                               25 contains symbolic info from crestr
c
      call tranvbini
      tvirtu = cpulft(0)
_IF(atmol)
c
c...  recl only valid for direct-access files
c
c     open(33,status='scratch',access='sequential',form='unformatted',
c    &      recl=300000)
      call fopen(33,'fort.33','priv','unformatted')
c     open(33,status='scratch',access='sequential',form='unformatted')
c
      call fopen(25,struct,'shared','unformatted')
_ENDIF
      write(iwr,12) struct
12    format(' && structures-file ',a,' &&',/)
c
c...  input
c
c...  vectors are now in q(1)
c...  at end of calculation ci-vector ends up at start of scratch space
c
         call izero(n35,noccurr,1)
c
         call wr15vb(ndum1,ndum2,ndum3,ndum4,ntstru,ndum6,ndum7,
     &               ndum8,ndum9,ndum10,ndum11,ndum12,'read')
         if (ovbci2) then
            nscf_old=nscf
            nscf=0
         endif
         if (nscf.le.0) then
c
c...        normal vb-calculation 
c
            kvlen = (max0(nbasis*nbasis,nbasbas*nbasbas))*2
            maxvec=kvlen/nbasis
            
            kv = igmem_alloc_inf(kvlen,'vbin.m','vb','kv',IGMEM_DEBUG)
            kci = igmem_alloc_inf(ntstru,'vbin.m','vb','kci',
     &                            IGMEM_DEBUG)
            call vbci(qq(kv),qq(kci),maxvec)
c
         else
c
c...        orbital optimalisation
c

            kvlen = (max0(nbasis*nbasis,nbasbas*nbasbas))*6
            maxvec = kvlen/nbasis
cjvl            write(iwr,1) kvlen,maxvec
1           format(/' *** vector array',i7,' for',i5,' vectors')
            kv = igmem_alloc_inf(kvlen,'vbin.m','vb','kv',IGMEM_DEBUG)
            kci = igmem_alloc_inf(ntstru,'vbin.m','vb','kci',
     +                            IGMEM_DEBUG)
c
c...        the array v in vbscf should be able to contain both
c...        the orbitals and the total virtual space
*
         if (oqcscf) then
*     
*           quadratic scf for orbital optimisation
*     
            call vbqcscf(qq(kv),maxvec,qq(kci),qq)
*    
         else
*
*           this is normal scf (super-ci or approx. newton-raphson)
*           for orbital optimisation
*
            call vbscf(qq(kv),maxvec,qq(kci))
*
         end if
*        
      end if
*        
      call gmem_free_inf(kci, 'vbin.m', 'vb', 'kci' )
      call gmem_free_inf(kv,  'vbin.m', 'vb', 'kv'  )
*
_IF(peigss)
         call cleargas()
_ENDIF
c
c...  that's it
c
c...  Show statistics of cases
c
         if (ovbci2) then
            nscf=nscf_old
         endif
      call showstat
c
c     call secsum
      call revind
c     call whtps
      call clredx
c
      total=tvirtu+t4indx0+t4indx+tmatre0+tmatre+tdavid+tlagr+tgtr+tfock
      nmatre0 = max0(nmatre0,1)
      n4indx0 = max0(n4indx0,1)
      nfock = max0(nfock,1)
c
      write(iwr,33)tmatre0,100.0d0*tmatre0/total,tmatre0/nmatre0,
     1             tmatre,100.0d0*tmatre/total,tmatre/nmatre,
     1             t4indx0,100.0d0*t4indx0/total,t4indx0/n4indx0,
     1             t4indx,100.0d0*t4indc/total,t4indx/n4indx,
     1             tfock,100*tfock/total,tfock/nfock,
     2             tvirtu,100.0d0*tvirtu/total,
     2             tdavid,100.0d0*tdavid/total,
     2             tlagr,100.0d0*tlagr/total,
     2             tgtr,100.0d0*tgtr/total,
     2             total
33    format(/
     &,' ============================================================',/
     &,'             timing analysis           sec    perc    per/it',/
     &,' matrix-elements psi0            :',1pe10.3,1x,0pf5.1,1pe11.4,/
     &,' matrix-elements orb-opt         :',1pe10.3,1x,0pf5.1,1pe11.4,/
     &,' integral transformation psi0    :',1pe10.3,1x,0pf5.1,1pe11.4,/
     &,' integral transformation orb-opt :',1pe10.3,1x,0pf5.1,1pe11.4,/
     &,' Fock matrix construction        :',1pe10.3,1x,0pf5.1,1pe11.4,/
     &,' virtual space construction      :',1pe10.3,1x,0pf5.1,/
     &,' diagonalisation                 :',1pe10.3,1x,0pf5.1,/
     &,' lagrangian and density matrices :',1pe10.3,1x,0pf5.1,/
     &,' density matrix transformation   :',1pe10.3,1x,0pf5.1,/
     &,' Total                           :',1pe10.3,1/
     &,' ===========================================================')
      write(iwr,34) minvar,maxvarp,minpert,maxpert,minfock,maxfock,
     &              minfockd,maxfockd
34     format('      Orbital optimisation methods'/
     &,' Var  : min',i9,'   max',i9,/' Pert : min',i9,'   max',i9,/
     &,' Fock : min',i9,'   max',i9,/' Fockd: min',i9,'   max',i9,/
     &,' ===========================================================')
c
      write(iwr,22) cpulft(1),charwall()
22    format(/,' leave program after ',f11.3,' CPU-seconds',a10,' wall')
c
      ofirst = .false.
c
      call flushbuffer
c
      if (irowhs.ne.1) then
         write(iwr,'(a)') 
     1   ' **** with incomplete H-matrix no use going on ****'
         call closda
_IF(parallel)
         call pg_end(0)
_ELSE
         stop    
_ENDIF
      end if
c
      return
      end

c***********************************************************************
      subroutine checkblk(words,block,tape,text)
c
      implicit REAL (a-h,o-z)
      common/checkblock/ icheckit
      dimension q(10000)
      character*(*) text
      integer words,block,tape
c

      if (icheckit .le. 0) return
      call secsum
      call rdedx(q,words,block,tape)

      return
      end

c***********************************************************************
      subroutine inivb
      implicit REAL  (a-h,o-z) , integer  (i-n)
c
c...  initialisation of vb-input
c
c=======================================================================
c
      common/bypass/ index4,ihmat,idavid
c
c     indicate (0=no/1=yes) to do resp. 4-index, matrix-el.,davidson
c
c=======================================================================
c
      common/title/ titel,chtime(4)
      character*80 titel
      character*8 chtime
c
c=======================================================================
c
INCLUDE(common/turtleparam)
INCLUDE(../m4/common/sizes)
      common/posit/iky(maxorb*(maxorb+1)/2)
c
c     large iky   for use in hsmat and getin2 from mclr
c
c=======================================================================
c
INCLUDE(common/hsinfo)
c
c     pass parameters from h-matrix generating part
c     nstruc   dimension of matrices
c     iprinv   print flag
c     ihfile,ihbloc   file and starting block for h and s-matrix
c     ipbloc  starting block for perturbation info
c
c=======================================================================
c
      common/davcom/eshift,thresj,thresd,maxcyc,maxdav,maxsel,iprind,
     *              alter,firsts
      common/davctx/mode
      character*4 mode
      logical alter,firsts
c
c     davidson variables
c
c=======================================================================
c
INCLUDE(common/brill)
INCLUDE(common/scftvb)
      common/davscf/eshsst,thrssj,thrssd,maxssc,maxssv,maxssl,iprssd,
     *              alssr,firsss
      common/davssx/modess
c
c=======================================================================
c
c     /cases/ is used for matrix-element type analysis (mtype)
      common /cases/ icase(2*ncases),times(ncases),token(3*ncases),
     &               lcases
      logical lcases
c
c======================================================================
c
      character*4 modess
      logical alssr,firsss
c
c     commons for orbital-optimisation control
c     nscf    occupied orbitals to optimise (0 => no opt)
c     iex    contains excitations
c            iex(1)   index of occupied orbital
c            iex(2..maxex)   indices of orbitals to mix in
c     nex     excitations for each orbital
c     nequi   equivalent orbital sets
c     iequi    starting adress per equivalence set in iex
c     ortscf   if .true. a.o.'s are orthogonalised onto occupieds
c
c     shiscf   scf-shift :(amount to substract from h11 in bi-diag)
c     criscf   convergence criterion (on orbital change)
c     evb    energy of vb-function
c     ebi    energy of bi-function
c     iprins   print-flag for scf
c     nitscf   iteration-count
c     maxscf   maximum  iterations
c     isecv    section for final vb-vectors (required input (?))
c
c     /davscf/ and /davssx/ are replacement of /davcom/,/davctx/ for bi
c     see subroutine davids for /davcom/ <=> /davscf/
c
c======================================================================
c
c...  /fast/ contains info on wether or not to use the symmetry due
c...         to equivalence relations in the matrix-element evaluation.
c...         (nosymm). nospin indicates wether spin-symetry may be used
c
      common /fast/ nospin,nosymm,nosymc
      logical nospin,nosymm,nosymc
INCLUDE(common/twice)
INCLUDE(common/infato)
c
c======================================================================
c...   ci and scf virtuals
INCLUDE(common/vbvirt)
c
c======================================================================
c
c...  /forbex/ will contain hand-given forbidden excitations.
c...           ment for eliminating redundancies in case psi(0)
c...           contains it's own brillouin state etc.
c...           also contains orthogonalisation info.
      common /forbex/ nforbi,iforbi(maxact*maxact,2),
     &                nionj,ionj(maxact,2),
     &                nsetor,isetor(maxact+1,maxact),
     &                nionres,ionrest(maxact)
c
c=======================================================================
c
      common/energi/enatom,nstate
c.....hand given energy of atoms, for print in 'natorb',nstate is the
c.....number of states the eigenvalues and -vectors will be printed
c.....(for the ci, or possibly the converged scf wavefunction)
c=======================================================================
c
      common /sdpert/ jstat,nrefst,dthres,e0
c.....control block for doubles selection in mrsdci calculation
c
c=======================================================================
c
      common/subinf/nwidth
c...  nwidt gives maximum width of matrix-element band written in
c...  HAMILT or BAMILT ; Is updated there
c=======================================================================
      logical opridet
      common/pridets/opridet
c...  opridet tells program to show matrix-elements in terms of
c...  determinants
c=======================================================================
INCLUDE(common/vbcri)
c...   common meant for all kind of (fixed) vb criteria
c...   critor .. criterium to determine orthogonality
c...   crilow .. criterium  for orthogonalization (e.g.jacoby in lowdin)
c...   cridep .. criterium for dependency (on c**2)
c...   symcri .. criterium for symmetry equivalence
c...   remcri .. criterium for including excitation on basis of 1-e ints
c...   crihyb .. max hybrid contam. allowed in matrix and (*100) vector
c...   cribri .. min. size of bcoeff to avoid annihilation
c...   crinor .. min. norm of a vector
c...   cribub .. a meaningful difference (in a bubble sort)
c...   criign .. size of matrix elements that really may be ignored
c...   cripop .. criterium for pipek/mizek localization
c...   normally critor = 1.0d-14, crilow = critor/10, cridep = critor*10
c...            symcri = critor*1000 remcri = 1.0d-14 crihyb = 1.0d-14
c...            cribri = 1.0d-8 crinor = critor*100  cribub = 1.0d-10
c...            criign = 1.0d-44 cripop 1.0d-10
c     REAL critor,crilow,cridep,symcri,remcri,crihyb,cribri,
c    1     crinor,cribub,criign,cripop
c     common/vbcri/ critor,crilow,cridep,symcri,remcri,crihyb,cribri,
c    1              crinor,cribub,criign,cripop
c=======================================================================
INCLUDE(../m4/common/dm)
c...  common block for gradients in GAMESS should be included 
c=======================================================================
INCLUDE(../m4/common/restar)
INCLUDE(../m4/common/restrj)
INCLUDE(../m4/common/restri)
INCLUDE(../m4/common/atmol3)
INCLUDE(common/vbparr)
INCLUDE(common/aivb)
c=======================================================================
c...   ed7 addressing for vb in kscra7vb
c     kbl7vb : current start of free space (must be updated)
c     k7struc : beginning of structures file (set in vbcrestr)
c     kvb7_bco,kvb7_vo,kvb7_vn: det-coef,old/new vectors 
c     kvb7_bcn                  set when first used in brigen and orbopt
c     kvb7_fcao,kvb7_fcmo     : fock-matrix on ao and mo bases
c     kvb7_do,kvb7_dn         : ao density matrix (old and new)
c     kvb1,kvb2,kvb3 : set when firsts used in diis
c=======================================================================
c
INCLUDE(common/hsmattype)
INCLUDE(common/vbproper)
      common /ncofac3/nul0,nul1,nul2,nul3,nul4,ndim3
c
c    /vbproper/
c
      orunprop=.false.
      oprop=.false.
      omag=.false.
      opol=.false.
      odebug=.false.
      onosym=.false.
      ovbci=.false.
      ocomplex=.false.
      ouncoupled=.false.
c
c
c,,,, default section numbers
c.... vectors 10
c.... natorb  11
c
      isecv = isect(310)
      if (mouta.le.0) mouta=isect(311)
c
c...  sync 1 : synchronous; check 1 : checking
c
      sync = 1
      check = 1
c
      irest = 0
c
c     /dm/
c
      ifmola = 1
      iflagr = 1
c
c
c     /bypass/
c
      index4 = 1
      ihmat = 1
      idavid = 1
c
c...  /posit/
c
      do 5 i=1,mxorbvb*(mxorbvb+1)/2
         iky(i) = i * (i-1) / 2
5     continue
c
c     /hsinfo/
c
      iprinv = 0
      ihfile = 6
      ihbloc = 1
      ipbloc = -1
c
c     /vbcri/
c
      critor = 1.0d-14
      crilow = 1.0d-15
      cridep = 1.0d-13
      crinor = 1.0d-12
      symcri = 1.0d-11
      remcri = 5.0d-15
c       remcri may not be too big (not d-10)
      crihyb = 1.0d-14
      cribri = 1.0d-8
      cribub = 1.0d-10
      criign = 1.0d-44
      cripop = 1.0d-10
c
c     /subinf/
c
      nwidth = 0
c
c     Initialise total number of vbscf runs 
c
      ntscf = 0 
      ipcmt = 0
c
c     /davcom/ , /davctx/
c
      eshift = 0.0d0
      thresj = 1.0d-14
      thresd = 5.0d-7
      maxcyc = 50
      maxdav = 30
      maxsel = 264
c     maxsel = 1
      iprind = 0
      alter = .false.
      firsts = .false.
      mode = 'emin'
c
c     /brill/ , /scf/
c
      nscf = 0
      do 10 j=1,maxact
         nex(j) = 0
         exact(j) = .true.
      do 10 i=1,maxex
10    iex(i,j) = 0
      nequi = 0
      ninter = 0
      nullity = 0
      ortscf = .true.
      reduce = .true.
      autscf = .false.
      schmi = .false.
      idequi = 0
      do i=1,maxact
         ieqs(i) = 0
         nexeq(i) = 0
      end do
c
      Optcri = 3
      shiscf = 0.0d0
      dampvb = 1.0d0
      fixb0 = -1.0d0
      criscf = 3.0d-5
      iprins = 0
      maxscf = 50
      call izero(maxex*maxact,ieqsig,1)
      call izero(maxact*(maxact+1),inteq,1)
      nredund = 0
      call izero(mxorbvb,iredund,1)
      maxdiis = 0
      percmax = 0.0d0
      perccur = 0.0d0
c
c     /davscf/ , /davssx/
c
      eshsst = 0.0d0
      thrssj = 1.0d-14
      thrssd = 1.0d-5
      maxssc = 50
      maxssv = 30
c     maxssl = 200
      maxssl = 1
      iprssd = 0
      alssr = .false.
      firsss = .true.
      modess = 'lock'

c
c     /cases/
c
      lcases = .false.
c
c      /atonam/ , /infato/
c
      call izero(maxato,nopa,1)
      call izero(maxato*maxopa,iopa,1)
      call izero(maxact*maxato,iacat,1)
      call izero(maxato,nacat,1)
      natom = 0
      mullik = .false.
      hybry = .false.
      clean = .true.
      hguess = .false.
      do 20 i=1,maxato
20    atoms(i) = '        '
c
c     /twice/
c
      call izero(maxact*5*maxact,inforb,1)
      call izero(maxact,ninfor,1)
      call izero(maxact,idoubly,1)
      ndoubly = 0
      call izero(maxact,isingly,1)
      nsingly = 0
      call izero(maxact*maxact,ioror,1)
      do 191 i=1,mxorbvb+maxact
191   kaart(i) = i
      super = .false.
      do i=1,maxact
         ofrzmo(i) = .false.
      end do
c
c     /fast/
      nospin = .false.
      nosymm = .false.
      nosymc = .false.
c
c     /vbvirt/
c     ci :
c     ortvirt causes the virtuals to be mutually orthogonal in a ci
c     calculation. if (provirt) the occupied orbitals are projected out
c     ivirt is the first virtual
c     scf :
c     canonicalise 2 : make orbitals diagonalise fock
c     localise : pipec localise - canonicalise=10
c     idempotent : make projection operator idempotent
c     aos : suggests to use ao's to mix in
c
      ortvirt = .false.
      provirt = .false.
      ivirt = 0
      canonicalise(1) = 0
      canonicalise(2) = 0
      idempotent = .false.
      aos = .false.
c
c...../forbex/ .
c
      nforbi = 0
      call izero(maxact*maxact*2,iforbi,1)
      nionj = 0
      call izero(maxact*2,ionj,1)
      nionres = 0
      call izero(maxact,ionrest,1)
      nsetor = 0
      call izero((maxact+1)*maxact,isetor,1)
c.....
c...../energi/
      enatom = 0.0d0
      nstate = 0
c
c...../sdpert/
c
      jstat = 0
      nrefst = 0
      dthres = 0.0d0
      e0 = 0.0d0
c
c..../pridets/
c
      opridet = .false.
c
c..../ncofac3/
c....keep track of nullity in third order cofactors for gradients
c....ment for debugging purposes
c
      nul0 = 0
      nul1 = 0
      nul2 = 0
      nul3 = 0
      nul4 = 0
      ndim3 = 0
c
c... restrj
c... set oclunp to true for unpack routine. In that case
c... the output array is cleared first
      oclunp = .true.
c
c...  /hsmattype/ (can be full or vbfunction)
c
      hsmatt='full      '
c
c...  p0core can be true or false
c
      p0core=.true.
c
      return
      end
      subroutine vbin(q,lword,qq)
      implicit REAL  (a-h,o-z) ,  integer   (i-n)
c
c...  *****************************************************************
c...           vb - 4index input routine
c...  *****************************************************************
c..
c..    vbin   returns vb (start) vectors  in q
c..    q must be big enough to hold 2 sets of vectors (for getqvb/vbin)
c
c...   communication with vb
c
      common/crestn/ struct
      character*44 struct
c
      common/bypass/ index4,ihmat,idavid
      common/title/ titel,chtime(4)
      character*80 titel
      character*8 chtime,ztest
      character*8 chyear
      character*1 ch1,ch2

c
c...  hsmat and bsmat
c
INCLUDE(common/hsinfo)
INCLUDE(common/vbcri)
c
c...   davidson
c
      common/davcom/eshift,thresj,thresd,maxcyc,maxdav,maxsel,iprind,
     *              alter,firsts
      common/davctx/mode
      character*4 mode
      logical alter,firsts 
_IF(parallel)
      logical oroot
_ENDIF
c
INCLUDE(common/c8_16vb)
c
c...   scf  (initialised in scfin)
c
INCLUDE(common/scftvb)
INCLUDE(common/turtleparam)
INCLUDE(common/vbparr)
INCLUDE(../m4/common/sizes)
INCLUDE(common/brill)
      common/davscf/eshsst,thrssj,thrssd,maxssc,maxssv,maxssl,iprssd,
     *              alssr,firsss
      common/davssx/modess
      character*4 modess
      logical alssr,firsss
c
c...   4-index common-blocks
c
INCLUDE(common/splice)
INCLUDE(common/basisvb)
INCLUDE(common/scra7vb)
INCLUDE(common/tractlt)
      common/junke/maxt,ires,ipass,nteff,npass1,npass2,
     *             nbuck,mloww,mhi,ntri,iacc
INCLUDE(../m4/common/mapper)
      common/discc/ ied(16)
      character*4 ied
c...  scratch  for mix-directive
      parameter (mxmix=20, mxorg=50)
      common /scrp/ cmix(mxmix,mxorg),imix(mxmix,mxorg),nopomx(mxorg)
c...
c
c...  cases for matrix type analysis
c
      common /cases/ icase(2*ncases),times(ncases),token(3*ncases),
     &               lcases
      logical lcases
c
INCLUDE(../m4/common/iofile)
c...
c...  atonam for names of atoms (for mulliken population analysis)
c...  nopa contains the # orbitals per atom, iopa contains them
c...  mullik indicates if the analysis is requested. iacat contains
c...  scf orbitals per atom, nacat the number of the, hybry indicates
c...  wether the automatic scf-procure should confine itself to
c...  atomic optomisation
c...
INCLUDE(common/infato)
INCLUDE(common/twice)
INCLUDE(common/vbvirt)
INCLUDE(common/dumpvbmo)
      integer ifchr,ilchr
      common/energi/enatom,nstate
INCLUDE(common/vbbonds)
c
c.....mrsdci
c     jstat control integer for doubles selection:
c
c     jstat < 0 : mrsdci without selection
c     jstat = 0 : no mrsdci (initialised in inivb)
c     jstat = 1 : mrsdci with doubles selection, phase 1
c     jstat = 2 : mrsdci with doubles selection, phase 2
      common /sdpert/ jstat,nrefst,dthres,e0
c
      logical opridet
      common/pridets/opridet
c        ... casa/casb
      common/actlen/lena(8),icfcor(100),ilifp(maxorb),len3(8),
     * ibext(8),lenext(8)
c
      logical dump,lschmi,llowdi
      common/comvbin/dump,lschmi,llowdi,nscset,nloset,ntemp,nmix
c
      logical new2el
      common/new2el/ new2el
c
      common/logvb/logcor,logact
      logical logcor,logact
c
c     MARCIN
      character*80 tmpfile,tmpdesc
c
      integer iprhs,irowhs
      common/hcpu_vb/iprhs,irowhs
      REAL prlev
      common/leading/prlev
c...
      dimension q(*),qq(*)
c...
      parameter (nkey=44 )
      character*4 idd(nkey),space,iby(3),iddp,test,ytest
c
c...   directive-keywords
c
      data idd/'pass','accu','acti','spli','onel',
     *         'vbve','basi','bypa','dipo','blks','curt','end ','fini',
     *         'titl','ipri','shif','crit','mix ','max ','case','sele',
     *         'mode','alte','scf ','mull','hybr','clea','nmos','schm',
     *         'lowd','moco','ener','eige','mrsd','dets','vbsc','para',
     *         'virt','hcpu','wpri','vbmo','new2','cres','chic'/
      data iby/'4ind','hmat','davi'/
      data iddp/'prin'/,space/' '/
c
c...  two scratch-files are used : 33 for vectors (length :
c...                               +/- 2*nbasis**2
c...                               25 contains symbolic info from crestr
c
      tvirtu = cpulft(1)
c
      logact = .false.
      logcor = .false.
c
c...  dump, to separate file, vb orbitals in proper format for molden
c...  first step, to making separate file, readable for molden with
c...  all the essential informations
c
      dumpvbmo = .false.
c
c...  new 2-el integrals ordering scheme by default off
c
      new2el = .false.
c
      iprhs = 0
      irowhs = 1
      nbonds = 0
c
c     call tidajt(chtime(2),chtime(3),chyear,chtime(1),chtime(4),isecs)
c     write(iwr,11) (chtime(ij),ij=1,4),isecs,lword
c     LOGO moved to vbfirst
_IF()
11    format(//,40x,50('_'),//
_IF(parallel)
     &,38x,'     t                      t       l              ',/
     &,38x,'     t                      t       l              ',/
     &,38x,'    tttt    u  u   rrrrr   tttt     l         eee  ',/
     &,38x,'     t      u  u    r   r   t       l        eeeee ',/
     &,38x,'     t   t  u  u    r       t   t   l   l    e     ',/
     &,38x,'      ttt    uu     r        ttt     lll      eee  '
     &,///
     &,9x, '  ----              ----       ',
     &     '     ----       ))))) --         ----   )))))      ',/
     &,9x, '       ----                   ',
     &     '----    _ _ _   { o o }  ----   _ _ _   { o o }    ',/
     &,9x, '                     ----     ',
     &     '      _(_)_(_)_  ( ^ ) ---    _(_)_(_)_  ( ^ )     ',/
     &,9x, ' ----         ----             ',
     &     '   _(_)_(_)_(_)_/ /        _(_)_(_)_(_)_/ /        ',/
     &,9x, '         ----              ----',
     &     ' _(_)_(_)_(_)_(_)/       _(_)_(_)_(_)_(_)/         ',/
     &,9x, '                               ',
     &     '(_)_(_)_(_)_(_)_(_)     (_)_(_)_(_)_(_)_(_)        ',/
     &,9x, ' ----     ----      ----       ',
     &     '  //  //   //  //         //  //   //  //           ',
_ELSE
     &,38x,'     t                      t       l              ',/
     &,38x,'     t                      t       l              ',/
     &,38x,'    tttt    u  u   rrrrr   tttt     l         eee  ',/
     &,38x,'     t      u  u    r   r   t       l        eeeee ',/
     &,38x,'     t   t  u  u    r       t   t   l   l    e     ',/
     &,38x,'      ttt    uu     r        ttt     lll      eee  '
     &,///
     &,    '  ----              ----       ',
     &     '     ----         ----           ----   )))))      ',/
     &,    '       ----                    ',
     &     '----       ----         ----   _ _ _   { o o }     ',/
     &,    '                     ----      ',
     &     '                   ----      _(_)_(_)_  ( ^ )      ',/
     &,    ' ----         ----             ',
     &     '      ----                 _(_)_(_)_(_)_/ /        ',/
     &,    '         ----              ----',
     &     '                         _(_)_(_)_(_)_(_)/         ',/
     &,    '                               ',
     &     '  ----                  (_)_(_)_(_)_(_)_(_)        ',/
     &,    ' ----     ----      ----       ',
     &     '           ----           //  //   //  //           ',
_ENDIF
     &          //40x,50('_'),
     &          //52x,'user        ',a8
     &          //52x,'date        ',a8
     &          //52x,'time        ',a8
     &          //52x,'job name    ',a8
     &          //52x,'job time    ',i6,' seconds',
     &          //52x,'main store',i8,' words')
c
c
c...  Create Structures (Formerly known as CRESTR)
c
_IF(atmol)
      call fopen(25,struct,'shared','unformatted')
_ENDIF
_ENDIF
      write(iwr,12) struct
12    format(' && structures-file ',a,' &&',/)
c
c...  input
c
      call flushbuffer
c
c
c.....some local initialisations for =schmidt=, =lowdin= and =hybrids=
      lschmi = .false.
      llowdi = .false.
      guesao = .false.
      nscset = 0
      nloset = 0
      prlev  = 0.001
c.....
      ncol = 0
      irest = 0
c*** TO BE DONE BY GAMESS
_IF(parallel)
      call setbfa(-1)
_ELSE
      call setbfa
_ENDIF
      call initra
      call inivb
      nmix = 0
      ntemp = 0
c.....ntemp possibly is the number of temporary mo's (see =nmos=)
c
c...   read   nbasis , iblk , dumpfile , section
c
c**   call input (already done in prep99)
c**   Now not done in crestr
c     call input
c 
c
c***  The following information should be supplied by GAMESS
c
c     call inpi(nbasis)
c     call inpi(iblk3)
c     call inpa4(test)
c     call inpi(ionsec)
c     iblk3 = max(iblk3,1)
      call cgamtovb(num3,iblk3)
      if (num3.le.0) num3 = 4
      if (ionsec.le.0) ionsec = 492
c
      write(iwr,6000) nbasis,ied(num3),iblk3,ionsec
6000  format(/' the number of basis functions is',i4//
     *' dumpfile on ',a4,'at block',i6,' 1e-integrals from section ',i3)
      if(nbasis.le.1)  call vberr(' dimension .le. 1 ')
c
      do i=1,mxorbvb
         ilifp(i)=(i-1)*nbasis
      end do
c
      lenbas=iky(nbasis+1)
      nbasbas = 0
c     call secini(iblk3,num3)
c     call secsum
c
c...   read directives
c
1     call input
      call inpa4(test)
      k=locatc(idd,nkey,test)
      if (k.eq.0) then
         call outrec
         call vberr('unrecognised directive')
      end if
c
      go to (5040,5050,5060,5070,5070,
     &       5090,5100,5110,5120,5130,5140,5150,5150,
     &       5170,5180,5190,5200,5210,5220,5230,5240,
     &       5250,5260,5270,5280,5290,5300,5310,5320,
     &       5330,5340,5350,5360,5370,5380,5270,5390,
     &       5400,5410,5420,5430,5440,5450,5460) , k

c...
c...   =pass= directive
c...
5040  call inpi(npass1)
      call inpi(npass2)
      if (npass1.lt.1.or.npass2.lt.1) call vberr('  passes too small')
      write(iwr,5041) ' *** multipass 4-index requested - passes ',
     1             npass1,npass2,' ***'
5041  format(a,2i6,a)
      goto 1
c...
c...   =accu= directive   // accuracy
c...
5050   call inpi(iacc)
       acc1 = 10.0d0**(-iacc)
       write(iwr,6050) acc1
6050   format(/' accuracy factor',e12.2)
       go to 1
c...
c...   =active=  directive
c...
5060  call inpa4(test)
      logact = test.eq.iddp
      call rdistr(mapiee,nsa,nsplice)
      if (nsa.eq.0) then
         call vberr(' no functions specified in active orbital list')
      else
         if (ncore_i.ne.0) then
            nsa = nsa - ncore_i
            do k=1,nsa
               if (mapiee(k).gt.ncore_i.and.k.le.ncore_i) 
     1            call caserr('core_i error in active')
               mapiee(k) = mapiee(k+ncore_i)
            end do
         end if
         write(iwr,6060) nsa
6060     format(/' the following',i5,'  mos included in active list'//
     *            5(' function  e/i label',4x))
         write(iwr,6061)(mapiee(k),k,k=1,nsa)
6061     format(5(i9,i11,4x))
         lenact=iky(nsa+1)
      end if
      go to 1
c ...
c ...=splice= or =onelec= directive
c ...
5070  call inpa4(test)
      if (ncore_i.ne.0) call caserr('splice *and* core')
      logcor = test .eq.iddp
      call rdistr(mapcie,ncore,nsplice)
      if (ncore.eq.0) then
         write(iwr,6070)
6070     format(/' no functions specified in frozen-core list'/)
      else
         write(iwr,6071)
6071     format(///' the following mos included in frozen core list'//
     *             5(' function  e/i label',4x))
         write(iwr,6061)(mapcie(k),k,k=1,ncore)
      end if
      go to 1
c...
c... =vbvectors= directive
c...
5090  call inpa4(test)
      isecvv = isecv
      if (test.eq.'manu') then
         write(iwr,5091)
5091     format(' vectors are supplied manually')
c.....   vectors given by hand
         call inpa4(test)
         dump = .false.
         if (test.eq.'dump') then
c.....      dump the vectors to be read on the dumpfile (only sensible
c.....      in a ci-type calculation
            dump = .true.
            call inpi(isecv)
            write(iwr,5099) isecv
5099        format(/,' write vectors to section ',i3)
         else
            call cinput(jrec,jump,-1)
            jrec = jrec - 1
         end if
         call inpi(ncoll)
         call inpa4(test)
         nnot = 0
         if (test.ne.'prin'.and.test.ne.'    ')  then
            call cinput(jrec,jump,-1)
            jrec = jrec - 1
            call inpi(nnot)
            call inpa4(test)
         end if
         call vecman(q(ncol*nbasis+1),nbasis,ncoll,nnot)
         ncol = ncol + ncoll
         if (test.eq.'prin') call prvc(q,ncol,nbasis,q,'o','l')
         if (dump) call putqvb(q,nbasis,ncol)
         isecv = isecvv
      else if (test.eq.'comb') then
c...     combine gamess section with vectors from input or other section
         call inpa4(test)
c...      
         maxcol = mxorbvb
         maxstr = 10
c
         if (3*maxcol*nbasis+maxcol*3+maxcol*maxstr*2.gt.lword)
     1      call caserr('not enough memory for veccomb')
         call veccomb(q(1),q(maxcol*nbasis+1),q(2*maxcol*nbasis+1),
     1                q(3*maxcol*nbasis+1),q(3*maxcol*nbasis+maxcol+1),
     2                q(3*maxcol*nbasis+maxcol*2+1),
     3                q(3*maxcol*nbasis+maxcol*3+1),
     4                q(3*maxcol*nbasis+maxcol*3+maxcol*maxstr+1),
     2                nbasis,ncol,maxcol,maxstr)
         if (test.eq.'prin') call prvc(q,ncol,nbasis,q,'o','l')
      else
          call cinput(jrec,jump,-1)
          jrec = jrec - 1
          call inpi(isec)
          call getqvb(q,nbasis,ncol,isec,'print')
          write(iwr,5092) isec
5092      format(' vectors restored from section',i3)
      end if
      irest = 1
      call inpa4(test)
      if (test .eq.iddp) call  prvc (q,ncol,nbasis,q,'o','l')
      go to 1
c...
c...  =basis=  isecb  ('print') : section that contains the basis to do the vb in 
c...
5100  call inpi(isecbas)
      call inpa4(test)
      isbass = isecbas
      call getqvb(q,nbasbas,ncolbas,isecbas,'print')
      write(iwr,5101) isecbas,nbasbas,nbasbas
5101  format(/1x,'     @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@',
     1       /1x,'     @@ VB is done in the basis of section',i6,' @@',
     2       /1x,'     @@ dimension',i6,' new dimension',i6,'      @@',
     3       /1x,'     @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@')
      if (nbasbas.ne.nbasis) call caserr('wrong basis dimension')
c...  fill scra7vb
      isecbas = 0
      k7bas =  kscra7vb('k7bas',nbasbas*ncolbas,'r','w')
      call wrt3(q,nbasbas*ncolbas,k7bas,num8)
      kfao = nbasbas*ncolbas+1
      kfmo = kfao + nbasbas*(nbasbas+1)/2
      ktemp = kfmo + ncolbas*(ncolbas+1)/2
      call get1e(q(kfao),dummy,'s',q(ktemp))
      call fmos(q(kfmo),q(kfao),q,q(ktemp),ncolbas,nbasbas,0.0d0)
      if (test.eq.'prin') then
         write(iwr,5102)
5102     format(/'           -----  new basis -----')
         call prvc (q,ncolbas,nbasbas,q,'o','o')
         write(iwr,5103)
5103     format(/'          -----  new metric -----')
         call prtri(q(kfmo),ncolbas)
         write(iwr,5104)
5104     format(1x)
      end if
      k7sbas =  kscra7vb('k7sbas',ncolbas*(ncolbas+1)/2,'r','w')
      call wrt3s(q(kfmo),ncolbas*(ncolbas+1)/2,num8)
      call get1e(q(kfao),dummy,'t',q(ktemp))
      call fmos(q(kfmo),q(kfao),q,q(ktemp),ncolbas,nbasbas,0.0d0)
      k7tbas =  kscra7vb('k7tbas',ncolbas*(ncolbas+1)/2,'r','w')
      call wrt3s(q(kfmo),ncolbas*(ncolbas+1)/2,num8)
      call get1e(q(kfao),dummy,'h',q(ktemp))
      call fmos(q(kfmo),q(kfao),q,q(ktemp),ncolbas,nbasbas,0.0d0)
      k7hbas =  kscra7vb('k7hbas',ncolbas*(ncolbas+1)/2,'r','w')
      call wrt3s(q(kfmo),ncolbas*(ncolbas+1)/2,num8)
      n7vbas = 0
c...  now in internal mode
      nbasis = ncolbas
      isecbas = isbass
      go to 1
c ...
c ... =bypass=   4index,hmatrix,davidson
c ...
5110  call inpa4(test)
      if (test.eq.iby(1)) index4 = 0
      if (test.eq.iby(2)) ihmat  = 0
      if (test.eq.iby(3)) idavid = 0
      call cinput(jrec,jump,0)
      if (jrec.lt.jump) go to 5110
      go to 1
c...
c... =dipole= directive
c...
5120  do 5121 i=1,3
5121  call inpf(dxyz(i))
      write(iwr,6120) dxyz
6120  format(/' dipole mixing factors - x=',f15.7,' y=',f15.7,
     *        ' z=',f15.7)
      if (isecbas.ne.0) call caserr('dipole should be before basis')
      go to 1
c...
c...  =blksize=  directive
c...
5130  call inpi(nsz)
      if (nsz.lt.1.or.nsz.gt.24) nsz=24
      write(iwr,6130) nsz
6130  format(/' sortfile blocking factor set to ',i3)
      go to 1
c...
c...   =curtail= directive  : limit first index of transformation
c...
5140  call inpa4(test)
      if (test.eq.'off') then
         ncurt = -1
      else if (test.eq.'on') then
         ncurt = nsa
      else
         call cinput(jrec,jump,-1)
         call inpi(ncurt)
         if (ncurt.eq.0) ncurt = nsa
         if (ncurt.eq.0) then
c...        pick up imax from crestr
_IF(atmol)
            rewind 25
            read(25) ndum1,ndum2,ndum3,ndum4,ndum5,ndum6,ndum7,
     &            ndum8,ndum9,imax
_ELSE
            call wr15vb(ndum1,ndum2,ndum3,ndum4,ndum5,ndum6,ndum7,
     &                  ndum8,ndum9,imax,ndum11,ndum12,'read')
_ENDIF
            ncurt = imax
         end if
         if (ncurt.eq.0) call vberr('cannot set ncurt')
      end if
      write(iwr,6140) ncurt
6140  format(//' @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@'
     *       /,' @   "curt"-option requested                    @'
     *       /,' @    first index curtailed to ',i6,'           @'
     *       /,' @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@')
      go to 1
c...
c...  =title= directive
c...
5170  read(ird,'(a)') titel
      write(iwr,6170) titel
6170  format(//,1x,132('|'),/,' |',t132,' |',/,' |',1x,a,t132,' |',/,
     &         ' |',t132,' |',/,1x,132('|'),/,1x)
      go to 1
c...
c...  =iprint=
c...
5180  call inpi(iiii)
      call cinput(jrec,jump,0)
      if (jrec.lt.jump) then
5185     call inpa4(test)
         if (test.eq.'davi') iprind = iiii
         if (test.eq.'tran') iprint = iiii
         if (test.eq.'hmat') iprinv = iiii
         if (test.eq.'scf')  iprins = iiii
         call cinput(jrec,jump,0)
         if (jrec.lt.jump) go to 5185
      else
         iprind = iiii
         iprint = iiii
         iprinv = iiii
         iprins = iiii
      end if
      write(iwr,6180) iprint,iprind,iprinv,iprins
6180  format(' iprint tran :',i5,' davidson :',i5,' hmat :',i5,' scf :',
     *       i5)
      go to 1
c...
c...  =shift=  davidson shift
c...
5190  call inpf(eshift)
      call outrec
      go to 1
c...
c...  =crit=  davidson (and/or) jacobi stop-criterion
c...          (and/or) orthogonality criteria
c...
5200  call inpf(xxx)
      call inpi(i)
      xxx = xxx * 10**(i*1.0d0)
      call cinput(jrec,jump,0)
      if (jrec.lt.jump) then
5201     call inpa4(test)
         if (test.eq.'davi') then
            thresd = xxx
            write(iwr,6200) thresd
         else if (test.eq.'jaco') then
            thresj = xxx
            write(iwr,6201) thresj
6201        format(/,' jacobi-criterion ',e12.5)
         else if (test.eq.'orth') then
            critor = xxx
            write(iwr,6202) critor
6202        format(/,' orthogonality-criteria +',e12.5)
            crilow = critor/10.0d0
            cridep = critor*10.0d0
            symcri = critor*100.0d0
            crinor = critor*100.0d0
         else if (test.eq.'lowd') then
            crilow = xxx
            write(iwr,6203) crilow
6203        format(/,' (lowdin) orthogonalisation-criterion ',e12.5)
         else if (test.eq.'depe') then
            cridep = xxx
            write(iwr,6204) cridep
6204        format(/,' dependency (on c**2) criterion ',e12.5)
         else if (test.eq.'symm') then
            symcri = xxx
            write(iwr,6205) symcri
6205        format(/,' symmetry (equiv) criterion ',e12.5)
         else if (test.eq.'exci') then
            remcri = xxx
            write(iwr,6206) remcri
6206        format(/,' excitation include criterion ',e12.5)
         else if (test.eq.'hybr') then
            crihyb = xxx
            write(iwr,6207) crihyb
6207        format(/,' hybrid check criterion ',e12.5)
         else if (test.eq.'bril') then
            cribri = xxx
            write(iwr,6208) cribri
6208        format(/,' minimum size of brillouin state ',e12.5)
         else if (test.eq.'igno') then
            criign = xxx
            write(iwr,6209) criign
6209        format(/,' minimum size of matrix elements ',e12.5)
         else if (test.eq.'pipe') then
            cripop = xxx
            write(iwr,6210) cripop
6210        format(/,' criterion for pipek-mizek localisation ',e12.5)
         end if
         call cinput(jrec,jump,0)
         if (jrec.lt.jump) go to 5201
      else
         thresd = xxx
         write(iwr,6200) thresd
6200     format(/,' davidson stop-criterion ',e12.5)
      end if
      go to 1
c...
c...  =mix =  mix input - orbitals "by hand"
c...
5210  continue
c
c...   if no orbitals have been given use unit matrix
c
      call inpa4(test)
      if (test.eq.'prin') logact = .true.
      if (irest.eq.0) then
         call extr(q,q,0,nbasis)
         irest = 1
         ncol = nbasis
         write(iwr,5211)
5211     format(/' -- no vectors are given => unit-matrix is used --')
      end if
5212  call input
      call inpa4(test)
      if (test.ne.'end') then
         call cinput(jrec,jump,-1)
         jrec = jrec - 1
         nmix = nmix + 1
         if (nmix.gt.mxorg) call vberr('too many orbitals in mix ')
         nopomx(nmix) = jump/2
         if (nopomx(nmix)*2.ne.jump) call vberr(' mix inconsistent')
         if (nopomx(nmix).gt.mxmix) call vberr(' to many mixes ')
         do 5213 i=1,nopomx(nmix)
            call inpi(imix(i,nmix))
            call inpf(cmix(i,nmix))
5213     continue
      else
         go to 1
      end if
      go to 5212
c...
c...  =max =  maximum  cycles in davidson   / max   expansion vectors
c...
5220  call inpi(nn)
      call cinput(jrec,jump,0)
      if (jrec.lt.jump) then
5221     call inpa4(test)
         if (test.eq.'davi') then
            maxcyc = nn
            write(iwr,6220) maxcyc
         else if (test.eq.'expa') then
            maxdav = min(250,nn)
            write(iwr,6221) maxdav
6221        format(' maximum  expansion-vectors in davidson ',i4)
         end if
         call cinput(jrec,jump,0)
         if (jrec.lt.jump) go to 5221
      else
         maxcyc = nn
         write(iwr,6220) maxcyc
6220     format(' maximum  davidson iterations ',i4)
      end if
      go to 1
c...
c...  =cases=  analyse the types of matrix-elements
c...
5230  lcases = .true.
      go to 1
c...
c...  =select=   states to select for starting vector
c...
5240  call inpi(maxsel)
      call inpa4(test)
      if (test.eq.'lowes') firsts = .false.
      maxsel = min(250,maxsel)
      if (firsts) then
         write(iwr,6240) ' first',maxsel
      else
         write(iwr,6240) 'lowest',maxsel
      end if
6240  format(' select ',a6,i5,' states for the start of davidson ')
      if (firsts) write(iwr,6241)
6241  format(' ** take the first not the lowest **')
      go to 1
c...
c...  =mode=   select  diagonalisation mode
c...
5250  call inpa4(test)
_IF(parallel)
      if (test.eq.'jacd') then
         write(iwr,*) '** Jacobi-Davidson not supported in parallel**'
         go to 1
      end if
_ENDIF
      mode = test
      if (mode.ne.'emin'.and.mode.ne.'vmin'.and.mode.ne.'lock'
     * .and.mode.ne.'jacd')
     *   call vberr(' wrong mode selected ')
      go to 1
c...
c...  =alter=  level-shift alteration in davidson
c...
5260  alter = .true.
      call inpa4(test)
      if (test.eq.'off') alter =.false.
      if (alter) write(iwr,6260)
6260  format(' level-shift alternation requested for davidson ')
      go to 1
c
c...  ==================================================================
c...
c...  =scf= or =vbscf= directives
c...
5270  if (nsa.eq.0) then
c
c...     get default nsa from crestr
c...     set nsa to imax / assumes crestr is called
c
_IF(atmol)
         rewind 25
         read(25) ndum1,ndum2,ndum3,ndum4,ndum5,ndum6,ndum7,
     &         ndum8,ndum9,imax
_ELSE
         call wr15vb(ndum1,ndum2,ndum3,ndum4,ndum5,ndum6,ndum7,
     &               ndum8,ndum9,imax,ndum11,ndum12,'read')
_ENDIF
         nsa = imax
         lenact=iky(nsa+1)
         do k=1,nsa
            mapiee(k) = k + ncore_i
         end do
         write(iwr,6060) nsa
         write(iwr,6061)(mapiee(k),k,k=1,nsa)
         lenact=iky(nsa+1)
      end if
      call scfin(q,lword,qq)
      go to 1
c...
c...  ==================================================================
c
c...
c...  =mulliken= atom definition in terms of orbital numbers
c...
5280  write(iwr,5282)
5282  format(/,' atom-definitions',/)
      mullik = .true.
5284  call input
      call inpa(ztest)
      if (ztest(1:4).eq.'end') go to 1
      natom = natom + 1
      atoms(natom) = ztest
      call inpa4(test)
      if (test.eq.'free'.or.test.eq.'frz') then
         ifrzat(natom) = 1
      else if (test.eq.'supe') then
         ifrzat(natom) = 2
      else
         ifrzat(natom) = 0
      end if
      call input
      call rdistr(iacat(1,natom),nacat(natom),maxact)
      if (ncore_i.gt.0) then
c...     eliminate core orbitals
5133     do i=1,nacat(natom)
            if (iacat(i,natom).le.ncore_i) then
               nacat(natom) = nacat(natom) - 1
               do j=i,nacat(natom)
                  iacat(j,natom) = iacat(j+1,natom)
               end do
               go to 5133
            end if
         end do
c...      renumber ..
         do i=1,nacat(natom)
            iacat(i,natom) = iacat(i,natom) - ncore_i
         end do
      end if
c...   freeze mo's
      if (ifrzat(natom).eq.1) then
         do i=1,nacat(natom)
            ofrzmo(iacat(i,natom)) = .true.
         end do
      end if
      call rdistr(iopa(1,natom),nopa(natom),maxopa)
      ii = natom
      test = '  '
      if (ifrzat(ii).eq.1) test = 'frz'
      if (ifrzat(ii).eq.2) test = 'sup'
      write(iwr,5286) atoms(ii),(iacat(ij,ii),ij=1,nacat(ii))
5286  format(1x,a8,1x,a3,' scf orbitals  ',(t28,20i4))
      write(iwr,5288) (iopa (ij,ii),ij=1,nopa(ii))
5288  format(1x,t10,' atomic orbitals  ',(t28,20i4))
      go to 5284
c.....
c.....=hybrids= try to make bonding hybrids according to rumer-bonds
c..... vor input or proximity
c.....
5290  hguess = .true.
      ichhyb = 0
      write(iwr,5291)
5291  format(/,' guess hybrids by maximising overlap')
5298  call inpa4(test)
      if (test.eq.'ao') then
         guesao = .true.
         write(iwr,*) ' per atom the aos will be used for this'
         go to 5298
      else if (test.eq.'atom') then
         ichhyb = 1
         write(iwr,*) ' bonded atom pairs are specified'
5297     call scann
         call inpa4(test)
         if (test.ne.'end') then
            call cinput(jrec,jump,-1)
            nbonds = nbonds + 1
            call rinpi(ibond(nbonds,1))
            call rinpi(ibond(nbonds,2))
            go to 5297
         end if
         go to 5298
      end if
      go to 1
c...
c...  =clear= clear mo's (i.e. remove coefficients of ao's from alien
c...          atoms. first atomic information has to be given through
c...          the mulliken option
c...
5300  if (.not.mullik) call vberr('please combine mulliken with clear')
5301  call input
      nmo = 0
5302  call inpa4(test)
      if (test.eq.'end') goto 1
      if (test.eq.'    ') call vberr('no atoms given in -clear-')
      ithere = 0
      do 5303 i=1,natom
5303  if (atoms(i).eq.test) ithere = 1
      if (ithere.eq.0) then
         nmo = nmo + 1
         goto 5302
      end if
      if (nmo.eq.0) call vberr('no mos given in -clear-')
      call cinput(jrec,jump,0)
      nat = jump - nmo
      do 5307 i=1,nmo
         call cinput(jrec,jump,i-1-jrec)
         call inpi(mocurr)
         do 5306 j=1,nat
            call cinput(jrec,jump,nmo+j-1-jrec)
            call inpa4(test)
            iat = 0
            do 5304 k=1,natom
5304        if (test.eq.atoms(k)) iat = k
            if (iat.eq.0) call vberr('unknown atom in -clear-')
            do 5305 k=1,nbasis
               ithere = locati(k,iopa(1,iat),nopa(iat))
               if (ithere.ne.0) q( (mocurr-1)*nbasis + k ) = 0.0
5305       continue
5306     continue
5307  continue
      goto 5301
c....
c.... =nmos= reduce the number of mo's defined in the active statement
c....        is usefull when mixing by hand. nsa2 has been set to zero
c....        at the beginning of this routine (used as flag)
c....
5310  call inpi(nsa2)
      if (nsa2.lt.1) call vberr('no active orbitals left over')
      write(iwr,5311) nsa2
5311  format(/,' the number of active orbitals will be reduced to ',
     &                                                     i3,/)
      ntemp = nsa - nsa2
      nsa = nsa2
      lenact = iky(nsa+1)
      goto 1
c.....
c.....=schmidt= orthogonalise vectors (is done just before leaving
c.....          this routine)
c.....
5320  lschmi = .true.
      nscset = nscset + 1
c.....use iex as scratch
      if (llowdi.and.nscset.gt.1)
     &   call vberr('too complicated, sorry')
      if (llowdi) write(iwr,5321)
5321  format(/,' === warning === schmidt will be done before lowdin ')
      call rdist1(iex(2,nscset),nnn)
      iex(1,nscset) = nnn
      write(iwr,5322) iex(1,nscset),
     &               (iex(ij,nscset),ij=2,iex(1,nscset)+1)
5322  format(/,' schmidt-orthogonalisation of the following ',i2,
     &         ' orbitals ',/,1x,(40i3),/,
     & ' N.B. numbering refers to state before "active" permutation')
      goto 1
c.....
c.....=lowdin= orthogonalise vectors (is done just before leaving
c.....         this routine)
c.....
5330  llowdi = .true.
      nloset = nloset + 1
c.....use iex as scratch
      if (lschmi) write(iwr,5331)
5331  format(/,' === warning === schmidt will be done before lowdin ')
      nbaal = nscset + nloset
      call rdist1(iex(2,nbaal),nnn)
      iex(1,nbaal) = nnn
      write(iwr,5332) iex(1,nbaal),(iex(ij,nbaal),ij=2,
     &                                     iex(1,nbaal)+1)
5332  format(/,' lowdin-orthogonalisation of the following ',i2,
     &         ' orbitals ',/,1x,(40i3),/,
     & ' N.B. numbering refers to state before "active" permutation')
      goto 1
c.....
c.....=mocopy= copy mo's just read (possibly permuted)
c.....
5340  if(ncol.eq.0)call vberr('please give vectors before -mocopy-')
      ipr = 0
5342  call inpa4(test)
      if (test.eq.'prin') ipr = 1
      call input
      call inpa4(test)
      if (test.eq.'end') then
         if (ipr.eq.1) call prvc(q,ncol,nbasis,q,'o','l')
         goto 1
      end if
      ncol = ncol + 1
      call cinput(jrec,jump,-1)
      jrec = jrec - 1
      call inpi(inew)
      call inpi(iold)
      if (iold.gt.ncol) call vberr('unknown old mo referenced')
c.....use iredund as scratch
      call rdistr(iredund,nao,mxorbvb)
      if(nao.ne.nbasis)call vberr('please specify all aos in mocopy')
      do 5341 i=1,nao
         if (iredund(i).ne.0) then
            ii = isign(1,iredund(i))
            iredund(i) = iabs(iredund(i))
            q((inew-1)*nbasis+i) = q((iold-1)*nbasis+iredund(i))*ii
         else
            q((inew-1)*nbasis+i) = 0.0d0
         end if
5341  continue
      goto 5342
c.....
c.....=energy= read sum of atomic energies (for print in "natorb")
c.....
5350  call inpf(enatom)
      write(iwr,5355) enatom
5355  format(/,' given energy of non-interacting atoms :',e21.15)
      goto 1
c.....
c.....=eigen= request print of ci eigenvalues and eigenvectors
 
5360  call inpa4(test)
      if (test.eq.'all') then
         nstate = 999999999
         write(iwr,5361)
5361     format(/,' print all eigenvalues and eigenvectors')
      else
         call cinput(jrec,jump,-1)
         jrec = jrec - 1
         call inpi(nstate)
         write(iwr,5362) nstate
5362     format(/,' print',i3,' lowest eigenvalues and eigenvectors')
      end if
      goto 1
c.....
c.....=mrsd= directive
c.....
5370  call inpa4(test)
      if (test.eq.'dsel') then
         call inpf(dthres)
         jstat =  1
      else
         jstat = -1
      end if
      go to 1
c....
c....=dets= requests output of matrix-elements between determinants
c....
5380  opridet = .true.
      go to 1
c....
c....=parallel=  sets parallel modes and checks
c....
5390  call inpa4(test)
      if (test.eq.'sync') then
         sync = 1
      else if (test.eq.'asyn') then
         sync = 0
      else if (test.eq.'chec') then
         check = 1
      else if (test.eq.'noch') then
         check = 0
      else
         if (sync.eq.1) write(iwr,5391) ' '
         if (sync.eq.0) write(iwr,5391) 'a'
         if (check.eq.1) write(iwr,5392) '  '
         if (check.eq.0) write(iwr,5392) 'no'
5391     format(' ** parallel mode : ',a1,'synchronous')
5392     format(' ** parallel mode : ',a2,' checksums')
         go to 1
      end if
      go to 5390
c....
c.... =virtual= ci virtuals directive
c....
5400  call caserr('ci virtual parameters may not yet be set')
      go to 1
c
c      =hcpu= print h/s matrix on the fly / allow arbitrary start
c             quick fix to make extremely timeconsuming vb restartable
c
5410  call inpi(irowhs)
      iprhs = 1
      irowhs = max(irowhs,1)
      write(iwr,5411) irowhs
5411  format(/' PRINT h/s matrix on the fly; start at group ',i6)
      go to 1
5420  call inpf(prlev)
      write(iwr,5421) prlev
5421  format(' weights above ',f13.10,' will be printed')
      go to 1
c
c...  =vbmo= dump vb mo orbitals to separate file
c
5430  continue
      dumpvbmo = .true.
      vbmounit = 20
      read(ird,'(a)') tmpdesc
      call stripblanks(ifchr,ilchr,tmpdesc)
      vbmodesc = tmpdesc(ifchr:ilchr)
      tmpfile = 'out.vbmo.'//vbmodesc
      vbmofile = tmpfile(1:len(tmpfile))
_IF(parallel)
      if ( oroot() ) then
        open(unit=vbmounit,status='new',file=vbmofile,form='formatted')
        goto 5431
      else
        goto 5432
      end if
_ENDIF
      open(unit=vbmounit,status='new',file=vbmofile,form='formatted')
      goto 5431
5432  write(iwr,'(a)') ' skipping node......'
5431  continue
      write(iwr,'(a)') ' '
      write(iwr,'(a)') ' VBMO: VB orbitals dumping requested '
      go to 1
c
c...  =new2= new ordering scheme of 2-el integrals
c
5440  continue
      new2el = .true.
      write(iwr,'(a)')
      write(iwr,'(a,a)') ' ** New 2-el integrals ordering scheme ',
     &                   'switched ON.'
      go to 1
c
c...  call crestr from vbin
c    
5450  call crestr(q(ncol*nbasis+1),lword-ncol*nbasis,q,qq)
      go to 1
c
c...  chicken ; tighten all criteria (On may undo this partially later)
c
5460  continue
      call inpa4(test)
      call inpa4(ytest)
      write(iwr,'(/,1x,a,a,a,a)') 
     1 ' **** chicken **** tighten criteria  ',test,'  ',ytest
c...   scftvb
      optcri = 3
      criscf = 3.0d-5
      maxscf = 50
c**   for prop criscf=min(criscf,1.0d-8)
c...   vbcri
      critor = 1.0d-14
      crilow = 1.0d-15
      cripop = 1.0d-10 
      cridep = 1.0d-13
      crinor = 1.0d-12
      symcri = 1.0d-11
c            symcri = critor*100.0d0
      remcri = 5.0d-15
      crihyb = 1.0d-14
      cribri = 1.0d-8
      cribub = 1.0d-10 
c...   davcom
      eshift = 0.0d0
      thresj = 1.0d-14
      thresd = 5.0d-7
      maxcyc = 50
      maxdav = 30
      maxsel = 264
c...   davscf
      eshsst = 0.0d0
      thrssj = 1.0d-14
      thrssd = 1.0d-5
      thrssd=min(thrssd,criscf*0.01d0)
      maxssc = 50
      maxssv = 30
      maxssl = 1
c
      if (test.ne.'off'.and.ytest.ne.'off') then
c
c***            chicken
c
c...   scftvb
         optcri = 3
         criscf = 1.0d-9
         maxscf = 250
c...   vbcri
         critor = 1.0d-40
         criign = 1.0d-44 
         crilow = 1.0d-40
         cripop = 1.0d-40
         cridep = 1.0d-40
         crinor = 1.0d-40
         symcri = 1.0d-40
         remcri = 5.0d-40
         crihyb = 1.0d-40
         cribri = 1.0d-40
         cribub = 1.0d-40
         criign = 1.0d-44
c...   davcom
         eshift = 0.0d0
         thresj = 1.0d-40
         thresd = 1.0d-14
         maxcyc = 200
         maxdav = 200
         maxsel = 264
c...   davscf
         eshsst = 0.0d0
         thrssj = 1.0d-40
         thrssd = 1.0d-30
         maxssc = 200
         maxssv = 200
         maxssl = 1
      end if
c
      if (test.eq.'prin'.or.ytest.eq.'prin') then
            write(iwr,'(a)') 
     1      ' ================= TURTLE criteria =================='
         if (cribri.eq.1.0d-40) then
            write(iwr,'(a)') ' TURTLE criteria, chicken is on'
         else
            write(iwr,'(a)') ' TURTLE criteria, chicken is off'
         end if
c...   scftvb
         write(iwr,'(a,i3)') ' optcri (div/oberlap/Brillouin) ',optcri
         write(iwr,'(a,1pe12.5)') ' criscf - stop for SCF ',criscf
         write(iwr,'(a,i10)') ' maxscf - max. # it. ', maxscf 
c...   vbcri
         write(iwr,'(a,1pe12.5)') 
     1    ' critor .. criterium to determine orthogonality ', critor 
         write(iwr,'(a,1pe12.5)')  
     1    ' criign .. size of matrix el. that may be ignored ',criign
         write(iwr,'(a,1pe12.5)') 
     1    ' crilow .. crit. for orthog. (e.g.jacoby in lowdin) ',crilow 
         write(iwr,'(a,1pe12.5)') 
     1    ' cripop .. crit for cizek-mezek localisation ',cripop 
         write(iwr,'(a,1pe12.5)')
     1    ' cridep .. criterium for dependency (on c**2) ',cridep
         write(iwr,'(a,1pe12.5)') 
     1    ' crinor .. min. norm of a vector ',crinor
         write(iwr,'(a,1pe12.5)') 
     1    ' symcri .. criterium for symmetry equivalence ',symcri
         write(iwr,'(a,1pe12.5)') 
     1    ' remcri .. crit. including excitation on 1-e ints ',remcri
         write(iwr,'(a,1pe12.5)') 
     1    ' crihyb .. max hybrid contam. in matrix and (*100) vector ',
     1    crihyb
         write(iwr,'(a,1pe12.5)') 
     1    ' cribri .. min. size of bcoeff to avoid annihilation ',cribri
         write(iwr,'(a,1pe12.5)') 
     1    ' cribub .. a meaningful difference (in a bubble sort) ',cribub
         write(iwr,'(a,1pe12.5)') 
     1    ' cripop .. criterium for pipek/mizek localization ',cripop
c...   davcom
         write(iwr,'(a,1pe12.5)') 
     1     ' eshift .. level shift for ci david ',eshift
         write(iwr,'(a,1pe12.5)') ' thresj .. jacobi crit ci ',thresj
         write(iwr,'(a,1pe12.5)') ' thresd .. davidson crit ci ',thresd
         write(iwr,'(a,i5)') 
     1    ' maxcyc .. max. # cycles davidson ci ',maxcyc
         write(iwr,'(a,i5)') 
     1    ' maxdav .. max. # vectors davidson ci ',maxdav
         write(iwr,'(a,i5)') 
     1    ' maxsel .. max. # vectors selected davidson ci ',maxsel
c...   davscf
         write(iwr,'(a,1pe12.5)') 
     1     ' eshsst .. level shift for scf david ',eshsst
         write(iwr,'(a,1pe12.5)') ' thrssj .. jacobi crit scf ',thrssj
         write(iwr,'(a,1pe12.5)') ' thrssd .. davidson crit scf ',thrssd
         write(iwr,'(a,i5)') 
     1    ' maxssc .. max. # cycles davidson scf ',maxssc
         write(iwr,'(a,i5)') 
     1    ' maxssv .. max. # vectors davidson scf ',maxssv
         write(iwr,'(a,i5)') 
     1    ' maxssl .. max. # vectors selected davidson scf ',maxssl
            write(iwr,'(a)') 
     1      ' ===================================================='
c
      end if
c
      go to 1
c...
c...  =end= or =finish= directive ... a vb is allowed
c...
5150  call inpa4(test)
      if (test.ne.'vb') call cinput(jrec,jump,-1)
      call inpi(i)
      if (i.ne.0) isecv = i
c
c...  get default nsa from crestr
c
      if (nsa.eq.0) then
c       set nsa to imax / assumes crestr is called
_IF(atmol)
         rewind 25
         read(25) ndum1,ndum2,ndum3,ndum4,ndum5,ndum6,ndum7,
     &         ndum8,ndum9,imax
_ELSE
         call wr15vb(ndum1,ndum2,ndum3,ndum4,ndum5,ndum6,ndum7,
     &               ndum8,ndum9,imax,ndum11,ndum12,'read')
_ENDIF
         nsa = imax
         lenact=iky(nsa+1)
         do k=1,nsa
            mapiee(k) = k + ncore_i
         end do
         write(iwr,6060) nsa
         write(iwr,6061)(mapiee(k),k,k=1,nsa)
      end if
      if (nsa.lt.1) call vberr(' no active orbitals')
c
c...  if no orbitals have been given use unit matrix
c
      if (irest.eq.0) then
         call extr(q,q,0,nbasis)
         ncol = nbasis
         write(iwr,604)
604      format(/,' ** no vectors are given ==> unit-matrix is used **')
      end if
c...
c...  check ncurt option
c...
      if (ncurt.eq.999999) then
         ncurt = nsa
         write(iwr,'(1x,a,i5,a)') '*** curtail set to ',ncurt,' ***'
      end if
      if (ncurt.gt.nsa) then
         write(iwr,605) ncurt,nsa
605      format(//' ***** ncurt (',i6,') greater than nsa (',i6,')')
         call vberr(' ncurt >  active orbitals ')
      end if
c
      write(iwr,6900) isecv
6900  format(/,' write vb-vectors to section ',i4)
      call putqvb(q,nbasis,ncol)
c
      nn = max(nbasis,nsa,ncol)
      if (nn.gt.mxorbvb) write(iwr,6950) nn,mxorbvb
6950  format(1x/1x,'**** max dimension found',i8,' max allowed'
     1       ,i8,' ****',/20x,'**** adapt turtleparam ****')
      if (nn.gt.255.and.n8_16.ne.16) write(iwr,6951) nn
6951  format(1x/1x,'**** max dimension ',i8,' >255',
     2       ' use a 255+ directive in CRESTR')
      if (nn.gt.mxorbvb.or.(nn.gt.255.and.n8_16.ne.16))
     1  call vberr(' dimension failure ')
c
c     call prvc(q,nsa+ncore+ntemp,nbasis,q,'o','l')
      return
      end
      integer function lchbas(v,q,nmo,action)
c
c...  do basis transformation
c...  action :
c...      external : switch to external basis (e.g. 4-index)
c...      internal : switch to internal basis
c...      off      : disable internal basis completely
c...                 also disable dipole field 
c...                 for SERVEC
c...      check    : check if another basis is in  effect
c...                 0  - no ; 1 - yes-active; -1 yes-inactive
c
      implicit REAL (a-h,o-z), integer (i-n)
      dimension v(*),q(*)
      character*(*) action
INCLUDE(common/basisvb)
INCLUDE(common/tractlt)
INCLUDE(common/splice)
INCLUDE(common/scra7vb)
INCLUDE(../m4/common/iofile)
      save n7nmo
c
      lchbas = 0
      if (isbass.ne.0) then
         lchbas = 1
         if (isecbas.eq.0) lchbas = -1
      end if
c
c...  OK if never internal basis but crash if in wrong state
c
      if (isbass.eq.0) return
c
      if (action.eq.'check') then
         return
      else if (action.eq.'external'.or.action.eq.'off') then
         if (isecbas.eq.0) call caserr('lchbas already external')
c...     save vectors in internal basis
         k7vbas = kscra7vb(k7vbas,nmo*nbasis,'r','w')
         n7vbas = nmo*nbasis
         call wrt3(v,nmo*nbasis,k7vbas,num8)
         nmo7bas = nmo
         nbas7bas = nbasis
c...     get vectors in external basis
         if (nbasis.ne.nbasbas) call caserr('lchbas check v-space')
         kqq = nbasbas*nbasis + 1
         call getqvb(q,nb,nc,isecbas,'noprint')
         if (nbasbas.ne.nb.or.ncolbas.ne.nc) call caserr('lchbas')
         if (ncolbas.ne.nbasis) call caserr('basis confusion')
c
         call mxma(q,1,nbasbas,v,1,nbasis,q(kqq),1,nbasbas,
     1             nbasbas,nbasis,nmo)
         call dcopy(nmo*nbasbas,q(kqq),1,v,1)
c
         nbasis = nbasbas
         isecbas = 0
         if (action.eq.'off') then
            isbass = 0
c...        switch off dipole directive as well
            do i=1,3
               dxyz(i) = 0.0d0
            end do
         end if
      else if (action.eq.'internal') then
         if (isecbas.ne.0) call caserr('lchbas already internal')
         if (n7vbas.eq.0) call caserr('we never were external')
         nbasis = ncolbas
         isecbas = isbass
         if (nmo.ne.nmo7bas) call caserr('nmo mismatch')
         if (nbasis.ne.nbas7bas) call caserr('nbas mismatch')
         k7vbas =  kscra7vb('k7vbas',nmo*nbasis,'r','r')
         call rdedx(v,nmo*nbasis,k7vbas,num8)
      else
         call caserr('wrong action')
      end if
c
      lenbas = nbasis*(nbasis+1)/2
c
      return
      end
   
      subroutine scfin(q,lword,qq)
c
      implicit REAL (a-h,o-z), integer (i-n)
c
      dimension q(lword),qq(*)
c
c...  input handling for scf
c...  called at scf directive / finished by end
c
INCLUDE(../m4/common/sizes)
INCLUDE(common/turtleparam)
INCLUDE(common/vbqc)
INCLUDE(common/brill)
c
c...  /forbex/ will contain hand-given forbidden excitations
c...  ment for eliminating redundancies in case psi(0)
c...  contains it's own brillouin state
c
      common /forbex/ nforbi,iforbi(maxact*maxact,2),
     &                nionj,ionj(maxact,2),
     &                nsetor,isetor(maxact+1,maxact),
     &                nionres,ionrest(maxact)
c
INCLUDE(common/twice)
INCLUDE(common/vbvirt)
c
      common/davcom/eshift,thresj,thresd,maxcyc,maxdav,maxsel,iprind,
     *              alter,firsts
      common/davctx/mode
      character*4 mode
      logical alter,firsts
c
      common /davscf/ eshsst,thrssj,thrssd,maxssc,maxssv,maxssl,iprssd,
     &                alssr,firsss
      logical alssr,firsss
c
      common /davssx/ modess
      character*4 modess
c
INCLUDE(common/tractlt)
INCLUDE(common/hsinfo)
INCLUDE(common/scftvb)
INCLUDE(common/vbcri)
INCLUDE(common/infato)
c
      common /fast/ nospin,nosymm,nosymc
      logical nospin,nosymm,nosymc
c
INCLUDE(common/vbdiist)
INCLUDE(common/vbpert)
c
      common /n2hybr/ n2hyb
      common /txema/ struc,ivec
      logical struc
c
INCLUDE(../m4/common/iofile)
      integer ifreez(maxact)
INCLUDE(common/splice)
INCLUDE(common/vbproper)
      integer itemp(mxorbvb)
INCLUDE(common/vbequiv)
      logical new2el
      common/new2el/ new2el
_IF(linux)
      external inpf
_ENDIF
c...
c...  =scf= directives
c...
      logical onot,ochatao,oflop,auto_hybrid
      parameter (nkey=38)
      character*4 idd(nkey),test,tes2,tes3,ytext
      character*8 ztest
      character*9 zfrom(4),zto(4)
      character*12 ztype(4)
      character*40 zstring
c
      data ztype/' super ci   ','perturbation',' ',' automatic  '/
      data zfrom/' doubly  ',' doubly  ','variably ','variably '/
      data zto  /'variably ','       un','variably ','      un '/
c
c...  directive-keywords
c
      data idd/'mix ','exci','crit','max ','mode','alte','ipri','shif',
     &         'sele','orth','remo','hybr','forc','nosy','nosp','supe',
     &         'eige','diis','forb','opti','evec','free','end','virt',
     &         'p0co','pcmt','prop','damp','leve','cano','fock','equi',
     &         'new2','curt','unit','qcsc','incs','scri'/
c
      auto_hybrid = .true.
      ovbdiis = .false.
      equiv = .false.
      equivir = .false.
      iwdiis = 0
      ntdiis = 9999999
      npert = 0
      ovbper = .false.
      fckp = .true.
      fockdiagonal = .false.
      ntype=0
      ifirsp = 99999
      idelp = -100
      ivec = 1
      oqcscf = .false.
      unitary = .false.
      scri = 0.0
c
1     call input
2     call inpa4(test)
      k = locatc(idd,nkey,test)
      if (k.le.0) call vberr('unrecognised directive')
c
      go to (5010,5010,5030,5040,5050,5060,5070,5080,
     &       5090,5100,5110,5120,5130,5140,5150,5160,
     &       5170,5190,5200,5210,5220,5230,5900,5240,
     &       5250,5260,5270,5280,5080,5290,5300,5310,
     &       5320,5330,5340,5350,5360,5370),k
c...
c... =excit=  or  =mix=
c...
5010  call inpa4(test)
      if (test.eq.'asis') then
         reduce = .false.
         write(iwr,6011)
6011     format(/' atomic orbitals will be mixed in exactly as',
     &              ' defined in the mix option',/,
     &              ' be shure not to mix in a doubly occupied part of',
     &              ' the orbital space !')
      else
         write(iwr,6012)
6012     format(/' atomic orbitals will be mixed in, while taking',
     &              ' the restriction in freedom, related to the pauli',
     &              ' principle, into account automatically')
      end if
5011  call input
      call inpa4(test)
      k = locatc(idd,nkey,test)
      call cinput(jrec,jump,-1)
      if (k.gt.0) go to 2
      if (test.eq.'equi') then
c...     equivalent orbital excitation patterns are defined
         call inpa4(test)
         call inpa4(test)
         if (test.eq.'inte') then
c....       equivalences of occupied orbitals are given (always last)
c....       use ninter as flag (negative value)
            ninter = -1
         else
            nequi = nequi + 1
            iequi(nequi) = nscf + 1
            if (test.eq.'exac') then
               exact(nequi) = .true.
            else
               exact(nequi) = .false.
            end if
         end if
         goto 5011
      end if
      if (ninter.eq.0) then
c...     no recognised directive / so should be orbital mix pattern
         nscf = nscf + 1
         if (nscf.gt.maxact)
     &                  call vberr('too many active orbitals in scf')
         call rdist1(iexw(1,nscf),nn)
         do 5013 i=1,nscf-1
            if (iexw(1,i).eq.iexw(1,nscf)) then
c...           extension of existing excitation pattern
               ione = 1
               do 5012 j=1,nn-1
                  ieqsig(nexw(i)+j,i) = isign(ione,iexw(j+1,nscf))
c....             sign-matrix for equivalence definitions
                  iexw(nexw(i)+j,i) = iabs(iexw(j+1,nscf))
5012           continue
               nscf = nscf - 1
               nexw(i) = nexw(i) + nn - 1
               if (nexw(i).gt.maxex)
     &                        call vberr('now too many ao s for mo')
               goto 5011
            end if
5013     continue
         if (nn.gt.maxex) call vberr(' too many ao s for 1 mo')
         nexw(nscf) = nn-1
         ione = 1
         do 5014 i=1,nn
            ieqsig(i,nscf) = isign(ione,iexw(i,nscf))
c....             sign-matrix for equivalence definitions
            iexw  (i,nscf) = iabs(iexw(i,nscf) )
5014     continue
      else
c....    used ninter as flag => reset it to zero
         if (ninter.lt.0) ninter = 0
         ninter = ninter + 1
         if (ninter.gt.nscf)
     &       call vberr(' equivalent definitions >   scf orbitals ?')
         call rdist1(inteq(1,ninter),nn)
         if (nn-1.gt.nscf)
     &       call vberr(' equivalent orbitals >   scf orbitals ?')
         ninteq(ninter) = nn - 1
      end if
      go to 5011
c...
c...  =crit=  for jacobi, davidson or scf (default scf)
c...
5030  call inpf(xxx)
      call inpi(i)
      xxx = xxx * 10**(i*1.0d0)
      call cinput(jrec,jump,0)
      if (jrec.lt.jump) then
5031     call inpa4(test)
         if (test.eq.'jaco') thrssj = xxx
         if (test.eq.'davi') thrssd = xxx
         if (test.eq.'scf')  criscf = xxx
         if (test.eq.'over') then
           optcri = 2
           criscf = xxx
           write(iwr,5333) criscf
5333       format(' use overlap-convergence criterion ',1pe7.1)
         end if
         if (test.eq.'bvec'.or.test.eq.'div ') then
           optcri = 1
           criscf = xxx
           write(iwr,5334) criscf
5334   format(' use largest brillouin-vector element criterion (div) ',
     1        1pe7.1)
         end if
         if (test.eq.'bmat') then
           optcri = 3
           criscf = xxx
           write(iwr,5335) criscf
5335   format(' use largest brillouin-matrix element criterion ',1pe7.1)
         end if
         call cinput(jrec,jump,0)
         if (jrec.lt.jump) go to 5031
      else
         criscf = xxx
      end if
      go to 1
c...
c...  =max=  for scf, davidson, expansion
c...
5040  call inpi(nn)
      call inpa4(test)
      if (test.eq.'scf '.or.test.eq.' ') then
         maxscf = nn
         write(iwr,5042) maxscf
5042     format(' scf : maximum # iterations ',i4)
      else if (test.eq.'davi') then
         maxssc = nn
         write(iwr,5043) maxssc
5043     format(' scf : maximum  davidson iterations ',i4)
      else if (test.eq.'expa') then
         maxssv = nn
         write(iwr,5044) maxssv
5044     format(' scf : maximum  davidson expasion vectors ',i4)
      endif
      call cinput(jrec,jump,0)
      if (jrec.lt.jump) go to 5040
c
      go to 1
c...
c...  =mode=
c...
5050  call inpa4(test)
_IF(parallel)
      if (test.eq.'jacd') then
         write(iwr,*) '** Jacobi-Davidson not supported in parallel**'
         go to 1
      end if
_ENDIF
      modess = test
      go to 1
c...
c...  =alter=
c...
5060  alssr = .true.
      go to 1
c...
c...  iprint  for scf and scf-davidson and the others
c...
5070  call inpi(iiii)
      call cinput(jrec,jump,0)
      if (jrec.lt.jump) then
5071     call inpa4(test)
         if (test.eq.'davi') iprssd = iiii
         if (test.eq.'scf')  iprins = iiii
         if (test.eq.'tran') iprint = iiii
         if (test.eq.'hmat') iprinv = iiii
         call cinput(jrec,jump,0)
         if (jrec.lt.jump) go to 5071
      else
         iprssd = iiii
         iprins = iiii
         iprint = iiii
         iprinv = iiii
      end if
      write(iwr,6070) iprssd,iprins,iprint,iprinv
6070  format(' scf iprint davidson  ',i5,'  scf  ',i5,' tran ',i5,
     1       ' hmat ',i5)
      go to 1
c...
c...   =shift=
c...
5080  call inpf(xx)
      call cinput(jrec,jump,0)
      if (jrec.lt.jump) then
5081     call inpa4(test)
         if (test.eq.'scf ') shiscf = xx
         if (test.eq.'davi') eshsst = xx
         call cinput(jrec,jump,0)
         if (jrec.lt.jump) go to 5081
      else
         shiscf = xx
      end if
      go to 1
c...
c...  select    (start for bi-diag) // defaults should be ok
c...
5090  call inpi(maxssl)
      call inpa4(test)
      if (test.eq.'firs') then
         firsss = .true.
      else if (test.eq.'lowe') then
         firsss = .false.
      end if
      go to 1
c...
c...   =ortho=
c...
5100  call inpa4(test)
      if (test.eq.'off') then
         ortscf = .false.
         write(iwr,5102)
5102     format(/' no orthogonalisations will be performed')
         goto 1
      end if
5103  call input
      call inpa4(test)
      if (test.eq.'end') then
         write(iwr,5111)
5111     format(/,' hand given orthogonolisation restrictions :',/)
         do 5104 i=1,nionj
5104     write(iwr,5144) ionj(i,1),ionj(i,2)
5144     format(1x,i3,' on',i3)
         do 5106 i=1,nionres
5106     write(iwr,5166) ionrest(i)
5166     format(1x,i3,' on the rest')
         do 5107 i=1,nsetor
5107     write(iwr,5178) (isetor(j,i),j=2,isetor(1,i)+1)
5178     format(' amongst one another',(t22,30i3))
         goto 1
      end if
      if (test.eq.'set') then
         nsetor = nsetor + 1
         call input
         call rdistr(isetor(2,nsetor),isetor(1,nsetor),maxact)
         if (nsetor.gt.maxact.or.isetor(1,nsetor).gt.maxact)
     &                 call vberr(' too many ortho-restrictions')
         goto 5103
      end if
      call cinput(jrec,jump,-1)
      call inpi(ii)
      call inpa4(test)
      call inpi(jj)
      if (test.eq.'on'.or.test.eq.'to') then
         nionj = nionj + 1
         if (nionj.gt.maxact) call vberr('too many ortho-restrictions')
         if (ii.le.0.or.jj.le.0) call vberr('wrong ortho input')
         ionj(nionj,1) = ii
         ionj(nionj,2) = jj
      else if (test.eq.'rest') then
         nionres = nionres + 1
         if(nionres.gt.maxact) call vberr('too many ortho-restrictions')
         if (ii.le.0) call vberr('wrong ortho input')
         ionrest(nionres) = ii
      else
         call vberr('wrong ortho input')
      end if
      goto 5103
c...
c...   =remove=  define smallness of one-electron int. as criterion
c...             in the automatic scf procedure
c...
5110  call inpf(xxx)
      call inpi(i)
      xxx = xxx * 10**(i*1.0d0)
      remcri = xxx
      write(iwr,6110) remcri
6110  format(' excitations corresponding to attraction integrals ',
     &       ' < ',e8.3,' will be removed')
      goto 1
c...
c...   =hybrid= define atoms for automatic atom-scf procedure.
c...
5120  call inpa4(test)
      if (test.eq.'clea') then
         clean = .true.
         go to 5120
      else if (test.eq.'dirt'.or.test.eq.'stra') then
         clean = .false.
         go to 5120
      else if (test.eq.'auto'.or.test.eq.'symm') then
         auto_hybrid = .true.
         go to 1
      else if (test.eq.'off') then
         auto_hybrid = .false.
         go to 1
      else if (test.eq.' ') then
c...   go on
      end if
c
      hybry = .true.
      auto_hybrid = .false.
c
      if (mullik) then
c...     atoms allready are defined
         write(iwr,5121)
5121     format(/,' the atom-definitions as given for the mulliken ',
     &             'population analysis are used for the hybrids too')
         goto 1
      end if
      natom = 0
5122  call input
      call inpa(ztest)
      if (ztest(1:4).eq.'gues') then
         write(iwr,5129)
5129     format(' spin-bonded hybrids will be guessed')
         hguess = .true.
         call input
         call inpa(ztest)
      end if
      if (ztest.eq.'end') then
c...  Within the hybrid option it is possible to leave out some  aos
c...  from all hybrid definitions. One could want to do this on purpose.
c...  Then an extra hybrid , carrying the name *rest* is created. 
c...  This hybrid has all remaining mo's and all unused aos.
         atoms(natom+1)="*rest*"
         ifrzat(natom+1) = 0
         nacat(natom+1)=0
         nopa(natom+1)=0
         do i=1,nbasis
            do j=1,natom
               do k=1,nopa(j)
                  if (iopa(k,j).eq.i) then
                     go to 6663
                  end if
               end do
            end do
            nopa(natom+1)=nopa(natom+1)+1
            iopa(nopa(natom+1),natom+1)=i
6663     continue
         end do
         if (nopa(natom+1).gt.0) then
            if (nsa.le.0) call caserr('no actives present for hybrid')
c....       add omitted mo's to this hybrid
            nacat(natom+1) = 0
            do 6664 i=1,nsa
              do j=1,natom
                do k=1,nacat(j)
                   if (iacat(k,j).eq.i) go to 6664
                end do
              end do
              nacat(natom+1) = nacat(natom+1) + 1
              iacat(nacat(natom+1),natom+1) = i
6664        continue
            natom=natom+1
         end if
c...
         write(iwr,5123)
5123     format(/,' hybrids will be made',/)
         if (clean) write(iwr,5127)
5127     format(' -- hybrid cleanness will be assured --')
         if (.not.clean) write(iwr,5128)
5128     format(' -- no cleaning of hybrids is performed --')
         do 5126 i=1,natom
            test = '   '
            if (ifrzat(i).eq.1) test = 'frz'
            if (ifrzat(i).eq.2) test = 'sup'
            write(iwr,5124) atoms(i),test,(iacat(ij,i),ij=1,nacat(i))
5124        format(1x,a8,1x,a3,' scf orbitals  ',(t28,20i4))
            write(iwr,5125) (iopa (ij,i),ij=1,nopa(i))
5125        format(1x,t10,' atomic orbitals  ',(t28,20i4))
5126     continue
         goto 1
      end if
      natom = natom + 1
      atoms(natom) = ztest
      call inpa4(test)
      if (test.eq.'free'.or.test.eq.'frz') then
         ifrzat(natom) = 1
      else if (test.eq.'supe') then
         ifrzat(natom) = 2
      else
         ifrzat(natom) = 0
      end if
      call input
      call rdistr(iacat(1,natom),nacat(natom),maxact)
      if (ncore_i.gt.0) then
c...     eliminate core orbitals
5133     do i=1,nacat(natom)
            if (iacat(i,natom).le.ncore_i) then
               nacat(natom) = nacat(natom) - 1
               do j=i,nacat(natom)
                  iacat(j,natom) = iacat(j+1,natom)
               end do
               go to 5133
            end if
         end do
c...     renumber ..
         do i=1,nacat(natom)
            iacat(i,natom) = iacat(i,natom) - ncore_i
         end do
      end if
      if (ifrzat(natom).eq.1) then
         do i=1,nacat(natom)
            ofrzmo(iacat(i,natom)) = .true.
         end do
      end if
      call input
c...   read atomic orbital definitions
      call inpa4(test)
      if (test.eq.'atom') then
c...     use atom defitions (and active mo's) to determine ao's
         kv = 1
         kvv = kv + ncol*nbasis
         kh = kvv + (ncol+nbasis)*nbasis
         kt = kh + nbasis*(nbasis+1)/2
         kmask = kt + (ncol+nbasis)*(ncol+nbasis+1)/2
         kihat = kmask + nbasis
         kscr = kihat + nbasis
         kend = kscr + nbasis
         if (kend.gt.lword) call caserr('core failure at atom_hybrid')
         call atom_hybrid(q(kv),q(kvv),q(kh),q(kt),q(kmask),q(kihat),
     1                    q(kscr))
      else 
c...     read ao's used in hybrid
         call cinput(jrec,jump,-1)
         call rdistr(iopa(1,natom),nopa(natom),maxopa)
      end if
c...  for large jobs where hybrids are defined as many-atom
c...  fragments the parameter maxopa might be too small
      if (nopa(natom).gt.maxopa) call vberr('maxopa too small')
      goto 5122
c.....
c.....=nosymm= tells not to use the symmetry.
c.....
5140  call inpa4(test)
      if (test.eq.'calc') then
         nosymc = .true.
         write(iwr,5145)
5145     format(/,' no spatial symmetry will be used in the matrix-'
     &        ,'element evaluation')
      else
         nosymm = .true.
         write(iwr,5147)
5147     format(/,' symmetry switched off')
      end if
      goto 1
c.....
c.....=nospin= tells not to use the symmetry due to spin
c.....in the matrix-element evaluation
c.....
5150  nospin = .true.
      write(iwr,5143)
5143  format(/' no spin-symmetry will be used in the matrix-element ',
     &         'evaluation')
      goto 1
c.....
c.....=super= refers to super-mcscf calculations, i.e. calculations in
c.....        which orbitals are optimised per csf (structure)
c.....        in that case no internal excitations are allowed for
c.....        as these cause numerical problems (they are replaced by
c.....        excitations to ao-s)
c.....
5160  super = .true.
      SUPER_hybrid = .false.
      super_cas = .false.
      super_act = .false.
      write(iwr,5165)
5165  format(/,' this is a super-mcscf calculation')
5161  call inpa4(test)
      if (test.eq.'hybr') then
         if ( natom .eq. 0 ) then
               call vberr(' unknown number of hybrids, use SUPER HYBRID
     &_AFTER_ HYBRIDS definition !!')
         end if
         super = .false.
         super_hybrid = .true.
         do i=1,natom
            igno_sh(i) = 1
         end do
         igno_sel = 2
c...     # of ao's to be ignored in generating excitations (from end)
c...     but if 1, we make sure it  appears in the orbital
5162     call inpa4(test)
         reduce = .false.
         if (test.eq.'igno'.or.test.eq.'excl') then
            do i=1,natom
               call inpi(igno_sh(i))
            end do
            go to 5162
         else if (test.eq.'last') then
            igno_sel = 1
            go to 5162
         else if (test.eq.'larg') then
            igno_sel = 2
            go to 5162
         else if (test.eq.'orth') then
            igno_sel = 99
            go to 5162
         else if (test.eq.'reduce') then
            reduce = .true.
            go to 5162
         end if
         nosymm = .true.
         nscf = nsa
         if (nsa.le.0) call vberr('active must be before super')
         if (igno_sel.eq.1) 
     1   write(iwr,5167) '(from end)',(igno_sh(i),i=1,natom)
         if (igno_sel.eq.2) 
     1   write(iwr,5167) '(largest )',(igno_sh(i),i=1,natom)
5167  format(' mixing is *completely* determined by hybrid directive',
     1     /,' ** symmetry switched off **',
     1     /,' ** AO ignored per set have to appear in MO (checked))',
     1     /,' ** AOs ',a10,' ignored :',20i4)
         if (igno_sel.eq.99)
     1   write(iwr,'(a)') ' generate seperate orth.compl. per hybrid'
         if (reduce) then
            write(iwr,'(a)') ' ** reduce is ON (BRRR) **'
         else
            write(iwr,'(a)') ' ** reduce is OFF (BRRR) **'
         end if
      else if (test.eq.'cas') then
         super_cas = .true.
         super = .false.
         write(iwr,5169) 
5169  format(' *no* active-active excitations allowed (CAS supposed)',/,
     1       ' Check the excitations; Optimisation may be incomplete ')
      else if (test.eq.'act'.or.test.eq.'acti') then
         super_act = .true.
         super = .false.
         write(iwr,5164) 
5164  format(' *no* excitations to virtuals allowed ',/,
     1       ' Check the excitations; Optimisation may be incomplete ')
      else
         goto 1
      end if
      go to 5161
c
c.... this one is free
c
5170  call vberr(' eigen is disabled')
c
c...  =diis= iwdiis ntdiis madiis midiis
c
5190  call inpi(iwdiis)
      call inpi(ntdiis)
      call inpi(madiis)
      call inpi(midiis)
      if (madiis.lt.3) madiis = 20
      if (midiis.lt.2) midiis = 3
      write(iwr,5195) iwdiis,ntdiis,madiis,midiis
5195  format(/' DIIS on if del < ',i2,' and',i2,' iterations',
     *       /' max # iterands',i4,' minimum # iterands',i4)
      go to 1
c.....
c.....=forbi= forbid (usually internal) excitations. ment for
c.....        elimination of redundancies in case psi(0) contains
c.....        it's own brillouin state. can be made automatic I
c.....        suppose.
5200  call inpi(ii)
      if (ii.ne.0) then
         nforbi = nforbi + 1
         iforbi(nforbi,1) = ii
         call inpi(ii)
         if (ii.eq.0) call vberr('inconsistent =forbi= input')
         iforbi(nforbi,2) = ii
         goto 5200
      else
         write(iwr,5201) (iforbi(i,1),iforbi(i,2),i=1,nforbi)
5201     format(/,' forbidden excitations :',(t25,15(i2,' =>',i2,2x)))
         goto 1
      end if
c.....
c..   =opti=  idelp  ifirsp  shiftp 
c.....
5210  call inpa4(test)
      call inpi(idelp)
      call inpi(ifirsp)
      call inpf(shiftp)
      if (test.eq.'bril') then
c...
c        =bril= optimisation per brillioun structure
c...
         write(iwr,*) 'scf optimisation specified per brillouin ',
     &   'structure (default super ci)'
         ntype=1
         do i=1,maxex
            itype(i)=0
         enddo
850      call input
         call inpa4(test)
         if (test.eq.'pert') then
            nn=2
            if ((ifirsp.le.0).or.(idelp.ge.1)) ovbper=.true.
         elseif (test.eq.'vari') then
            nn=1
         elseif (test.eq.'end ') then
            goto 860
         else
            call vberr('invalid =opti= input')
         endif
         call rdist1(itemp,itot)
         do i=1,itot
c... Ground state is brillouinstate 1. We want to be quasi zero based (like
c... the brillouin print) so we add 1 to every index
            if (itype(itemp(i)+1).eq.0) then
               itype(itemp(i)+1)=nn
            else
               write(iwr,*) 'Error: doubly specified orbital'
               call vberr('invalid =opti= input')
            endif
         enddo
         if (itot.le.30) then
            write(iwr,858) ztype(nn),(itemp(i),i=1,itot)
858         format(' using ',a12,' for structures ',30i3)
         else
            write(iwr,859) ztype(nn)
859         format(' using ',a12,' for more than 30 orbitals!')
         endif
         goto 850
860      do i=1,maxex
           if (itype(i).eq.0) itype(i)=1
         enddo
      elseif (test.eq.'orbi') then
c...
c        =orbi= optimisation per orbital
c...
         write(iwr,*) 'scf optimisation specified per orbital ',
     &   '(default super ci)'
         ntype=2
         do i=1,mxorbvb
            itype(i)=0
         enddo
870      call input
         call inpa4(test)
         if (test.eq.'pert') then
            nn=2
            if ((ifirsp.le.0).or.(idelp.ge.1)) ovbper=.true.
         elseif (test.eq.'vari') then
            nn=1
         elseif (test.eq.'end ') then
            goto 880
         else
            call vberr('invalid =opti= input')
         endif
         call rdist1(itemp,itot)
         do i=1,itot
            if (itype(itemp(i)).eq.0) then
               itype(itemp(i)-ncore_i)=nn
            else
               write(iwr,*) 'Error: doubly specified orbital'
               call vberr('invalid =opti= input')
            endif
         enddo
         if (itot.le.30) then
            write(iwr,878) ztype(nn),(itemp(i),i=1,itot)
878         format(' using ',a12,' for orbitals ',30i3)
         else
            write(iwr,879) ztype(nn)
879         format(' using ',a12,' for more than 30 orbitals!')
         endif
         goto 870
880      do i=1,mxorbvb
           if (itype(i).eq.0) itype(i)=1
         enddo
      elseif (test.eq.'kind') then
c...
c        =kind= optimisation per kind brillioun structure
c...
         ntype=3
         do i=1,4 
            itype(i)=1
         enddo
c...
890      call input
         call inpa4(test)
         call inpa4(tes2)
         call inpa4(tes3)
         select case(test)
           case('vari') 
              nn = 1
           case('pert') 
              nn = 2
           case('fock')
              call caserr('fock disabled')
           case('auto')
              nn = 4
           case('end')
           write(iwr,894)
           do i =1,4
             if (itype(i).eq.4) write(iwr,895) zfrom(i),zto(i),
     &                            ztype(itype(i)), crivarpert(i)
             if (itype(i).ne.4) write(iwr,896) zfrom(i),zto(i),
     &                                         ztype(itype(i))
           end do
894        format(1x,/,' === excitation control by kind  ===')
895        format('  from ',a10,'occupied to ',a10,
     &           'occupied treated by ',a12,' crit ',1pe10.3)
896        format('  from ',a10,'occupied to ',a10,
     &           'occupied treated by ',a12)
           go to 1
           case default
              call vberr('optimisation type unknown')
         end select
         if (nn.eq.4) then
            call inpf(crivarper)
            if (crivarper.eq.0.0d0) crivarper = 0.0001
         end if
         select case(tes2)
           case('doc')
             if (tes3.eq.'voc'.or.tes3.eq.'all') then
                itype(1) = nn
                if (nn.eq.4) crivarpert(1) = crivarper
             end if
             if (tes3.eq.'uoc'.or.tes3.eq.'all') then
                itype(2) = nn
                if (nn.eq.4) crivarpert(2) = crivarper
             end if
           case('voc')
             if (tes3.eq.'voc'.or.tes3.eq.'all') then
                itype(3) = nn
                if (nn.eq.4) crivarpert(3) = crivarper
             end if
             if (tes3.eq.'uoc'.or.tes3.eq.'all') then
                itype(4) = nn
                if (nn.eq.4) crivarpert(4) = crivarper
             end if
           case('all')
             if (tes3.eq.'voc'.or.tes3.eq.'all') then
                itype(1) = nn
                if (nn.eq.4) crivarpert(1) = crivarper
             end if
             if (tes3.eq.'voc'.or.tes3.eq.'all') then
                itype(3) = nn
                if (nn.eq.4) crivarpert(3) = crivarper
             end if
             if (tes3.eq.'uoc'.or.tes3.eq.'all') then
                itype(2) = nn
                if (nn.eq.4) crivarpert(2) = crivarper
             end if
             if (tes3.eq.'uoc'.or.tes3.eq.'all') then
                itype(4) = nn
                if (nn.eq.4) crivarpert(4) = crivarper
             end if
           case default
             call vberr('optimisation set unknown')
         end select
c
         if ((ifirsp.le.0).or.(idelp.ge.1)) ovbper=.true.
c
         goto 890
c
      else
         call vberr('invalid =opti= input')
      endif
c     else
c        call cinput(jrec,jump,-1)
c     end if
      if (ifirsp.le.0.or.idelp.ge.1) ovbper = .true.
      if (.not.ovbper) then
         write(iwr,5222) idelp,ifirsp
5222     format(/,' perturbation approach for virtuals after del ',i2,
     *           ' and iteration ',i2)
      else
         write(iwr,5233)
5233     format(/,' perturbation theory applied from the start ')
      end if
      if (shiftp.ne.0.0d0) write(iwr,5244) shiftp
5244  format(' level shift for PT ',f9.5)
      go to 1
c...
c...  =evec=
c...
c
c...  For optimising the orbitals of excited states
c...  By default, turtle optimises the orbitals for
c...  the lowest eigenvector in the basis of struc-
c...  tures (evec = 1). For optimising the orbitals
c...  of excited states choose the appropriate vec-
c...  tor with the evec directive. The order of the
c...  vectors may change during optimisation, which
c...  should be checked by the user.
c
5220  call inpi(ivec)
      goto 1
c
c... =freez= request to freeze mo's
c
5230  call rdistr(ifreez,nfreez,maxact)
      write(iwr,5231) nfreez,(ifreez(i),i=1,nfreez)
5231  format(1x,/' **freeze** ',i5,' mos ',(3x,t25,10i5),1x)
      do i=1,nfreez
         ofrzmo(ifreez(i)) = .true.
      end do
      go to 1
c
c... =virtual= directive to determine scf-virtuals
c              cf common/vbvirt
c
5240  onot = .false.
5241  call inpa4(test)
      if (test(1:2).eq.'no') then
        onot = .true.
      else if (test.eq.'cano') then
        if (onot) then
           canonicalise(2) = 0
        else
           canonicalise(2) = 1
        end if
        onot = .false.
      else if (test.eq.'fock') then
        if (onot) then
           canonicalise(2) = 0
        else
           canonicalise(2) = 2
           fckp = .true.
        end if
        onot = .false.
      else if (test.eq.'loca') then
        if (.not.onot) canonicalise(2) = 10
        if (canonicalise(2).eq.10) then
          nit_loc = 100
          call inpa4(test)
          if (test.eq.'nit'.or.test.eq.'iter') then
             call inpi(nit_loc)
          else
             call cinput(jrec,jump,-1)
          end if
        end if
        onot = .false.
      else if (test.eq.'idem') then
        idempotent = .not.onot
        onot = .false.
      else if (test.eq.'aos') then
        aos = .not.onot
        onot = .false.
        if (aos) canonicalise(2) = 0
      else
        write(iwr,'(a)') ' == Virtual orbital SCF options =='
        if (canonicalise(2).eq.1) write(iwr,'(a)') 
     1    ' Virtuals in SCF are canonicalised over h'
        if (canonicalise(2).eq.10) write(iwr,'(a)') 
     1         ' Virtuals in SCF are localised '
        if (idempotent) write(iwr,'(a)') 
     1    ' Projection operators are made idempotent'
        if (aos) write(iwr,'(a)') ' We try to use aos as virtuals '
        go to 1
      end if
      go to 5241
c...
c...  =p0co= - Psi(0) core seting can be "on" or "off". When Psi(0) is
c...           built the doubly occupied orbitals can be treated as core
c...           ("on") or as active ("off"). When p0core is not specified
c...           the default will be "on".
c...
5250  call inpa4(test)
      if (test.eq.'on') then
         p0core=.true.
      else if(test.eq.'off') then
         p0core=.false.
      else
         call vberr('invalid =p0co= input')
      endif
      write(iwr,5251) test
5251  format(/,' doubly treated as core for psi-0 p0-core-option --',a3)
      goto 1
c...
c...  =pcmt= - PCM tweak option. When the VB wave functions for the
c...           gasphase and the solution differ immensely, the VBPCM
c...           model is unable to converge to the right solution.
c...           With this option the structure number that will be
c...           important in the solution can be chosen to be important
c...           in the gasphase curve as well (tweak), which is the
c...           first VBSCF run in a PCM job.
c...
c...  example1: pcmt 2 5
c...           Run a 5 iterations VBSCF run with the emphasis on struc-
c...           ture 2. If the 5 is not specified VBSCF will run until
c...           the convergence criterion is met.
c...
c...  example2: pcmt 1000 1
c...           Run a 1 iteration VBSCF run (one iteration does not
c...           change the orbitals) leaving the vb-vector intact (caused
c...           by specifying 1000 or more).
c...
5260  call inpi(ipcmt)
      call inpi(maxscf_s)
      if (maxscf_s.eq.0) then
         maxscf_s = maxscf
      endif
      goto 1
5270  oprop=.true.
      criscf=min(criscf,1.0d-8)
      thrssd=min(thrssd,criscf*0.01d0)
5271  call inpa4(test)
      if (test.eq.'pola') then
          opol=.true.
      else if (test.eq.'magn') then
          omag=.true.
      else if (test.eq.'debu') then
          odebug=.true.
      else if (test.eq.'nosy') then
          onosym=.true.
      else if (test.eq.'igno') then
          oignore=.true.
      else if (test.eq.'vbci') then
          ovbci=.true.
      else if (test.eq.'unco') then
          ouncoupled=.true.
      else if (test.eq.'vbc2') then
          ovbci2=.true.
          call inpi(ndubbel)
      else if (test.eq.'end') then
          goto 1
      else
          call vberr('invalid =prop= input')
      endif
      goto 5271
c...
c...  =damp=
c...
5280  call inpf(xxxx)
      call inpa4(test)
      if (test.eq.'fixb'.or.test.eq.'b0') then
         fixb0 = xxxx
      else
         dampvb = xxxx
      end if
      goto 1
c
c...  canonicalise
c
5290  call inpa4(test)
      if (test.eq.'off') then
          canonicalise(1) = 0
          canonicalise(2) = 0
      else if (test.eq.'both') then
          call inpa4(test)
          iii = 0
          if (test.eq.'h') iii = 1
          if (test.eq.'fock') iii = 2
          if (test.eq.'loca') iii = 10
          canonicalise(1) = iii
          canonicalise(2) = iii
      else if (test.eq.'doub'.or.test.eq.'doc'.or.test.eq.'occ') then
          call inpa4(test)
          iii = 0
          if (test.eq.'h') iii = 1
          if (test.eq.'fock') iii = 2
          if (test.eq.'loca') iii = 10
          canonicalise(1) = iii
      else if (test.eq.'virt'.or.test.eq.'uoc') then
          call inpa4(test)
          iii = 0
          if (test.eq.'h') iii = 1
          if (test.eq.'fock') iii = 2
          if (test.eq.'loca') iii = 10
          canonicalise(2) = iii
      else if (test.eq.' ') then
         write(iwr,'(1x/a)') ' Canonicalisation requested'
         if (canonicalise(1).eq.0) write(iwr,'(a)') 
     1      ' Doubles over nothing'
         if (canonicalise(1).eq.1) write(iwr,'(a)') ' Doubles over h'
         if (canonicalise(1).eq.2) write(iwr,'(a)') ' Doubles over fock'
         if (canonicalise(1).eq.10) write(iwr,'(a)')' Doubles localised'
         if (canonicalise(2).eq.0) write(iwr,'(a)') 
     1      ' Virtuals over nothing'
         if (canonicalise(2).eq.1) write(iwr,'(a)') ' Virtuals over h'
         if (canonicalise(2).eq.2) write(iwr,'(a)') 
     1      ' Virtuals over fock'
         if (canonicalise(2).eq.10) write(iwr,'(a)') 
     1      ' Virtuals localised'
         go to 1
      end if
      if (canonicalise(1).eq.10.or.canonicalise(2).eq.10) then
         nit_loc = 100
         call inpa4(test)
         if (test.eq.'nit'.or.test.eq.'iter') then
            call inpi(nit_loc)
         else
            call cinput(jrec,jump,-1)
         end if
      end if
      go to 5290
c
c.... =fock=
c
5300  call inpa4(test)
      if (test.eq.'on'.or.test.eq.'yes') then
         fckp = .true.
      else if (test.eq.'off'.or.test.eq.'no') then
         fckp = .false.
      else if (test.eq.'diag') then
         fockdiagonal = .true.
      else if (test.eq.'nodi') then
         fockdiagonal = .false.
      else if (test.eq.' ') then
         if (fckp) write(iwr,'(2a)') 'Fock matrix used for Brillouin',
     1                              ' where appropriate'
         if (.not.fckp) write(iwr,'(2a)') 'no Fock matrix used for ',
     1                                    'Brillouin'
         if (fockdiagonal) write(iwr,'(2a)') 'Fock approximation for ',
     1                                       'diagonal with fock'
         if (.not.fockdiagonal) write(iwr,'(2a)') 'no Fock approxi',
     1                                            'mation for diagonal'
         go to 1
      end if
      go to 5300
c...
c...  =Equivalence force or  primitive= 
c...
5310  call inpa4(test)
      if (test.eq.'forc') then
         go to 5130
c...     we may go there directly
      else if (test.eq.'prim'.or.test.eq.'virt'.or.test.eq.' ') then
         if (test.eq.'virt'.or.test.eq.' ') call cinput(jrec,jump,-1)
         go to 5319
      else
         call caserr('wrong equiv keyword')
      end if
c
c....
c.... =equiv force= force group of scf-orbital2 into equivalence eg.px/y in o2
c....         should be used if orbitals are not equivalent
c....         in the csf's (structures), but are not so only because
c....         of the orbital-model. this option is ment for the auto-scf
c....         mode. the number refers to the group
c....         **note** equivalence orbitals are next to each other
c
5130  call rdistr(iequi,nequi,maxact)
      idequi = 1
5131  call inpa4(test)
      if (test.eq.' ') go to 5135
      if (test.eq.'asis') idequi = 1
c...   select excitations as done before
      if (test.eq.'h') idequi = 2
c...   select excitations on maximum h-matrix element
      if (test.eq.'n') idequi = 3
c...   do not select excitations, i.e. take the next one
      if (test.eq.'inte') then
c...   specify internal excitations - for an equivalence group
c...   pairs for each equivalece are given; finish with 'end'
         isss = 0
         do i=1,maxact
            isss = isss + nexeq(i)
         end do
         call rinpi(ieq)
         ieqs(ieq) = isss+1
5132     call inpa4(test)
         if (test.ne.'end') then
            call cinput(jrec,jump,-1)
c...        read equivalent excitations
            if (isss+2*iequi(ieq).gt.maxact*maxact) 
     &          call vberr('iexeq overflow')
            do i=isss+1,isss+2*iequi(ieq)
               call rinpi(iexeq(i))
            end do
            isss = isss + 2*iequi(ieq)
            nexeq(ieq) = nexeq(ieq) + 2*iequi(ieq)
            go to 5132
         end if
      end if
      go to 5131
5135  nn = 0
      do iforce=1,nequi
         exact(iforce) = .false.
         nn = nn + iequi(iforce)
      end do
      do is=nequi+1,nequi+nsa-nn 
         iequi(is) = 1
      end do
      nequi = nequi + nsa-nn
c      write(iwr,*) 
      write(iwr,'(/,a,20I3)') ' forced equivalence sets: ',
     1(iequi(is),is=1,nequi)
      if (idequi.eq.2) write(iwr,*) ' virtuals ordered by h-interaction'
      if (idequi.eq.3) write(iwr,*) ' virtuals not selected '
      do i=1,maxact
         if (nexeq(i).gt.0) then
            write(iwr,'(a,i3)') ' hand given excitations for class ',i
            isss = ieqs(i)-1
            nn = nexeq(i) / (2*iequi(i))
            do j=1,nn
               write(iwr,5136) (iexeq(k+isss),k=1,2*iequi(i))
5136           format(1x,(8(i4,' - ',i4,3x),1x),1x)
               isss = isss + 2*iequi(i)
            end do
         end if
      end do
      goto 1
c
c...  Equivalence primitive - may be united with the symmtry equiv
c...  only one set is allowed right now / only benzeen
c...  input is like hybrid ....
c
5319  equiv = .true. 
c
      call inpa4(test)
      if (test.eq.'virt') equivir = .true.
      neqmo = 0
5311  call input
      call inpa(ztest)
      if (ztest.eq.'end') then
         write(iwr,5312)
5312     format(/,' ***** primitive equivalences requested *****')
         if (equivir) write(iwr,5316)
5316     format(  '       generate equivalent virtuals  ')
         do i=1,neqmo
            write(iwr,5313) eqvmo(i)
5313        format(1x,' scf orbital  ',(t16,20i4))
            write(iwr,5315) (eqvao(ij,i),ij=1,neqao)
5315        format(1x,t10,' atomic orbitals  ',(t28,20i4))
         end do
         goto 1
      end if
      call cinput(jrec,jump,-1)
      neqmo = neqmo + 1
      if (neqmo.gt.meqmo) call caserr('equiv overflow')
      call inpi(eqvmo(neqmo))
c...   read atomic orbital definitions
      call input
      call inpa4(ztest)
      if (ztest.eq.'atom') then
c...     use atom defitions (and active mo's) to determine ao's
         kv = 1
         kvv = kv + ncol*nbasis
         kh = kvv + (ncol+nbasis)*nbasis
         kt = kh + nbasis*(nbasis+1)/2
         kmask = kt + (ncol+nbasis)*(ncol+nbasis+1)/2
         kihat = kmask + nbasis
         kscr = kihat + nbasis
         kend = kscr + nbasis
         if (kend.gt.lword) call caserr('core failure at atom_hybrid')
         call atom_hybridn(q(kv),q(kvv),q(kh),q(kt),q(kmask),q(kihat),
     1                     eqvmo(neqmo),1,eqvao(1,neqmo),nneqao,meqao,
     2                     q(kscr))
      else 
c...     read ao's used in equivlence
         call cinput(jrec,jump,-1)
         call rdistr(eqvao(1,neqmo),nneqao,meqao)
      end if
      if (neqmo.eq.1) then
         neqao = nneqao
         if (neqao.gt.meqao) call caserr('to many equiv ao')
      else
         if (nneqao.ne.neqao) call caserr('nonequal equiv #')
      end if
      goto 5311
c...
c...  =new2= new ordering scheme of 2-el integrals
c
5320  new2el = .true.
      write(iwr,'(/a,a)') ' ** New 2-el integrals ordering scheme ',
     &                   'switched ON.'
      go to 1
c...
c...   =curtail= directive  : limit first index of transformation
c...
5330  call inpa4(test)
      if (test.eq.'off') then
         ncurt = -1
      else if (test.eq.'on') then
         ncurt = nsa
      else
         call cinput(jrec,jump,-1)
         call inpi(ncurt)
         if (ncurt.eq.0) ncurt = nsa
      end if
      write(iwr,5331) ncurt
5331  format(//' @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@'
     *       /,' @   "curt"-option requested                    @'
     *       /,' @    first index curtailed to ',i6,'           @'
     *       /,' @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@')
      go to 1
c...
c...  unitary asks for an orthogonal mcscf
c...
5340  unitary = .true.
      write(iwr,*) '*** an orthogonal mcscf is requested ***'
      go to 1
*
*     = qcsc = directive (quadratic scf procdeure requested)
*
5350  oqcscf = .true.
      ocong = .true.
      ospci = .false.
      oinve = .false.
      ograd = .false.
      ocorr = .false.
      oauto = .false.
      oiter = .false.
      oqcof = .false.
      nitsci = 0
5351  call input
      call inpa4(test)
      select case(test)
         case('spci')
            ospci = .true.
            oqcof = .true.
            sci2nr = 0.01d0
            call inpa4(tes2)
            if (tes2.eq.'auto') then
               call inpf(sci2nr)
               call inpa4(tes3)
               if (tes3.eq.'grad') then
                  ograd = .true.
               else if (tes3.eq.'corr') then
                  ocorr = .true.
               else if (tes3.eq.'auto'.or.tes3.eq.' ') then
                  oauto = .true.
               end if
            else if (tes2.eq.'iter') then
               oiter = .true.
               call inpi(nitsci) 
            end if
         case('augm')
            call inpf(scalhs)
         case('inve')
            oinve = .true.
            ocong = .false.
         case('cong')
            ocong = .true.
            oinve = .false.
         case('end')
*
            write(iwr,5352)
5352        format(1x,89('*'))
            write(iwr,5353)
            if (ospci) then
               if (oiter) write(iwr,5354) nitsci
               if (ograd) write(iwr,5355) sci2nr
               if (ocorr) write(iwr,5356) sci2nr
               if (oauto) write(iwr,5357) sci2nr
            end if
5353        format(/,9x,
     +          'quadratically convergent scf procedure requested for '
     +          'orbital optimisation')
5354        format(15x,'use super-ci in the first ', i3,' iterations')
5355        format(15x,'use super-ci until the orbital gradient drops',
     +             ' to ',f10.5)
5356        format(15x,'use super-ci until the orbital correction'
     +             ' drops to ',f10.5)
5357        format(5x,'use super-ci until either the orbital gradient'
     +             ' or the correction vector drops to ',f10.5)
            write(iwr,5358)
5358        format(/,1x,89('*'))
*
            go to 1
*
         case default
            call vberr('unknown vbscf directive')
         end select
*
      go to 5351
*
c...
c...  =incs=
c...
5360  incsort = .true.
      write(iwr,*) ' *** incore sorting is on ***' 
       go to 1
c...
c...  scri
c...
5370  scri = 0.7
      call inpa(ytext)
5371  if (ytext.eq.'print') then
         oprs = .true.
      else
         jrec = jrec-1
         call inpf(scri)
         if (scri.eq.0.0d0) scri = 0.7
         if (ytext.ne.' ') go to 5371
      endif
      print *,'scri set to',scri
      go to 1
c...
c...  =end=
c...
5900  call inpa4(test)
      	if (test.ne.'scf'.and.test.ne.'vbsc') then
           call cinput(jrec,jump,-1)
        end if
      call inpi(isecv)
c
      if (nscf.eq.0) then
         write(iwr,6902) remcri
6902     format(/,' the automatic scf-procedure is invoked, criterion  '
     &           ,e8.3)
c        autscf = .true.
c...      autscf does not seem to have a function left **check**
         nscf = nsa
_IF(atmol)
         rewind 25
         read(25) ndum1,ndum2,ndum3,ndum4,ndum5,ndum6,ndum7,
     &         ndum8,ndum9,imax
_ELSE
         call wr15vb(ndum1,ndum2,ndum3,ndum4,ndum5,ndum6,ndum7,
     &               ndum8,ndum9,imax,ndum11,ndum12,'read')
_ENDIF
         if (imax.ne.nscf) then
            write(iwr,6903)
6903        format(/,' the number of active orbitals must be equal to ',
     &                'the highest occupied orbital number')
            call vberr(' try again')
         end if
      end if
c
      if (auto_hybrid) then
         call toetoec_hybrid(q,nsa,nbasis,qq)
         write(iwr,'(1x,/,1x,a,/,1x,a)') 
     1         'Automatic hybrid defnition',
     2         'switched off by hybrid off'
         if (clean) write(iwr,5127)
         if (.not.clean) write(iwr,5128)
         do i=1,natom
            write(iwr,5124) atoms(i),'   ',(iacat(ij,i),ij=1,nacat(i))
            write(iwr,5125) (iopa (ij,i),ij=1,nopa(i))
         end do
         hybry = .true.
      end if
c
c...  check hybrids and set atomao
c
      n2hyb = 0
      call izero(nbasis,atomao,1)
      do i=1,natom
        do j=1,nacat(i)
          if (iacat(j,i).gt.nsa) then
             write(iwr,6904) i,iacat(j,i),nsa
6904         format(' ** for atom',i3,' mo',i3,' gt nactiv (',i3,')')
             call vberr(' error in hybrid specification ')
          end if
        end do
      end do
      oflop = ochatao(flop,flop,'set')
c 
c check compatibility convergence criteria
c
      thrssd=min(thrssd,criscf*0.01d0)
      thresd=min(thresd,criscf*0.01d0)
      if (optcri.eq.2) then
        zstring = 'wavefunction overlap (over)'
      else if (optcri.eq.1) then
        zstring = 'Brillouin Interaction vector (bvec/div)'
      else if (optcri.eq.3) then
        zstring = 'Brillouin Theorem (bmat)'
      else 
         zstring = '*unknown*'
      end if
      write(iwr,7000)criscf,zstring,thresd,thrssd
 7000 format(1x,/,1x,88(1h*),/,1x,12x,' ** VB convergence criteria **'/
     1' ** SCF stop criterion - criscf      ',1pe8.2,' on ',a40,/
     2' ** Davidson Criterion (CI) - thresd ',1pe8.2,/
     3' ** Davidson criterion (BI) - thrssd ',1pe8.2,/,1x,88(1h*))
      if (shiscf.ne.0.0d0) write(iwr,7001) shiscf
      if (dampvb.ne.1.0d0) write(iwr,7002) dampvb
      if (fixb0.ge.0.0d0) write(iwr,7003) fixb0
 7001 format(' ** Level Shift for SCF              ',f8.3)
 7002 format(' ** Damping for SCF                  ',f8.3)
 7003 format(' ** b0 fixed in SCF to               ',f8.3) 
      return
      end
      subroutine toetoec_hybrid(q,nmdim,ndim,qq)
      implicit none
c
c...  calls actual worker
c
INCLUDE(../m4/common/gmemdata)
INCLUDE(../m4/common/gmempara)
      REAL qq(*),q(*)
      integer nmdim,ndim,kmoao,kim1,nipw,igmem_alloc_inf,kcc,krest
c
      kmoao = igmem_alloc_inf((ndim*nmdim+1)/nipw(),'vbin.m',
     &                        'toetoe_hybrid','kmoao',IGMEM_DEBUG)
      kim1 = igmem_alloc_inf((nmdim+1)/nipw(),'vbin.m',
     &                        'toetoe_hybrid','kim1',IGMEM_DEBUG)
      kcc = igmem_alloc_inf(ndim*(ndim+nmdim),'vbin.m',
     &                      'toetoec_hybrid','kcc',IGMEM_DEBUG)
      krest = igmem_alloc_inf((ndim+1)/nipw(),'vbin.m',
     &                        'toetoec_hybrid','krest',IGMEM_DEBUG)
c
      call toetoe_hybrid(q,nmdim,ndim,qq(kcc),qq(kmoao),qq(kim1),
     1                   qq(krest),qq)
c
      call gmem_free_inf(krest,'vbin.m','toetoec_hybrid','krest')
      call gmem_free_inf(kcc,'vbin.m','toetoec_hybrid','kcc')
      call gmem_free_inf(kim1,'vbin.m','toetoec_hybrid','kim1')
      call gmem_free_inf(kmoao,'vbin.m','toetoec_hybrid','kmoao')
c
      return
      end
      subroutine toetoe_hybrid(c,nmdim,ndim,cc,imoao,im1,irest,qq)
      implicit none
c
c.... set hybrids automatically to avoid symmetry contamination
c
INCLUDE(../m4/common/gmemdata)
INCLUDE(../m4/common/gmempara)
INCLUDE(common/turtleparam)
INCLUDE(common/vbcri)
INCLUDE(common/infato)
INCLUDE(common/splice)
      integer ndim,nmdim,nnmdim,imo,ii,i,nnn,ih,jmo
      integer kc,khao,khmo,kv,igmem_alloc_inf,locati,kel
      REAL dum
      REAL qq(*),c(ndim,*),cc(ndim,ndim+nmdim)
      integer imoao(ndim,nmdim),im1(nmdim),irest(ndim)
      nnmdim = ndim + nmdim
c
c...   set occupied + ao's
c
      call vclr(cc,1,ndim*nnmdim)
      call fmove(c(1,ncore+1),cc,nmdim*ndim)
      do imo=1,ndim
         cc(imo,imo+nmdim) = 1.0d0
      end do
c
c...   get the 1-electron matrix and transform to mo-basis
c
      khao = igmem_alloc_inf(ndim*(ndim+1)/2,'vbin.m',
     &                       'toetoe_hybrid','khao',IGMEM_DEBUG)
      call get1e(qq(khao),dum,'h',qq(khao)) 
      khmo = igmem_alloc_inf(nnmdim*(nnmdim+1)/2,'vbin.m',
     &                       'toetoe_hybrid','khmo',IGMEM_DEBUG)
      kv = igmem_alloc_inf(ndim,'vbin.m','toetoe_hybrid','kv',
     &                     IGMEM_DEBUG)
      call fmos(qq(khmo),qq(khao),cc,qq(kv),nnmdim,ndim,crilow)
c
c...   determine seperate groups (hybrids)
c
      call izero(nmdim*ndim,imoao,1)
      do imo=1,nmdim
         do jmo=nmdim+1,nnmdim
            kel = khmo-1 + jmo*(jmo-1)/2 + imo
            if (dabs(qq(kel)).gt.symcri) imoao(jmo-nmdim,imo) = 1
         end do
      end do
c
      call gmem_free_inf(kv,'vbin.m','toetoe_hybrid','kv')
      call gmem_free_inf(khmo,'vbin.m','toetoe_hybrid','khmo')
      call gmem_free_inf(khao,'vbin.m','toetoe_hybrid','khao')
c
c...  merge same groups
c
      do imo=1,nmdim
         do jmo=1,nmdim
            ii = 0
            do i=1,ndim
               ii = ii + imoao(i,imo)*imoao(i,jmo)
            end do
            if (ii.ne.0) then
               do i=1,ndim
                  if  (imoao(i,imo).ne.0) imoao(i,jmo) = 1
                  if  (imoao(i,jmo).ne.0) imoao(i,imo) = 1
               end do
            end if
         end do
      end do
c
c...  determine hybrids
c
      do imo=1,nmdim
         im1(imo) = locati(1,imoao(1,imo),ndim)
      end do
      call izero(ndim,irest,1)
c
      natom = 0
      nnn = 0
1     natom = natom + 1
      ih = 0
      if (natom.gt.maxato) call caserr('many hybrids in toetoe_hybrid')
      write(atoms(natom),'(a6,i2)') 'hybrid',natom
      ifrzat(natom) = 0
      nacat(natom) = 0
      do imo=1,nmdim
         if (im1(imo).lt.0) cycle
         if (ih.eq.0) then
            ih = imo
            nopa(natom) = 0
            do jmo=im1(imo),ndim
               if (imoao(jmo,imo).eq.1) then
                  nopa(natom) = nopa(natom) + 1
                  iopa(nopa(natom),natom) = jmo
                  irest(jmo) = 1
               end if
            end do
         end if
         if (iabs(im1(ih)).ne.iabs(im1(imo))) cycle
         nacat(natom) = nacat(natom) + 1
         iacat(nacat(natom),natom) = imo
         im1(imo) = -iabs(im1(imo))
         nnn = nnn + 1
      end do
      if (nnn.lt.nmdim) go to 1
      if (nnn.gt.nmdim) call caserr(' too many in toetoe_hybrid')
c
c...  gather rest
c
      imo = 0
      do jmo=1,ndim
         imo = imo + irest(jmo)
      end do
      imo = ndim - imo
      if (imo.gt.0) then
         natom = natom + 1 
         if (natom.gt.maxato) call caserr('many hybrids in toetoe_hyb')
         write(atoms(natom),'(a6,i2)') '*rest*',natom
         nacat(natom) = 0
         nopa(natom) = 0
         do jmo=1,ndim
            if (irest(jmo).eq.0) then
               nopa(natom) = nopa(natom) + 1
               iopa(nopa(natom),natom) = jmo
             end if
         end do
         if (nopa(natom).ne.imo) call caserr('oeps toetoe_hybrid')
      end if
c
      return
      end
      subroutine scfina(v,ci,s,q,qq)
c
      implicit REAL (a-h,o-z), integer(i-n)
c
c...  final print for scf
c
c...  if converged, return scfconv = .true.
c...  v    all active orbitals (no core)
c...  ci   ci vector
c...  s    mo overlap
c...  q    scratch
c
      dimension v(*),ci(*),s(*),q(*),qq(*)
c 
c...  for dynamic memory allocation debugging purpose, contains IGMEM_ vars.
INCLUDE(../m4/common/gmemdata)
INCLUDE(../m4/common/gmempara)
c
INCLUDE(common/dumpvbmo)
INCLUDE(common/turtleparam)
INCLUDE(../m4/common/sizes)
INCLUDE(../m4/common/common)
      common /energi/ enatom,nstate
c
INCLUDE(common/tractlt)
INCLUDE(common/ffile)
INCLUDE(common/hsinfo)
c
c...  nstruc (hsinfo) is now accidently dimension of bi(*)
c
INCLUDE(common/scftvb)
INCLUDE(common/splice)
INCLUDE(common/twice)
INCLUDE(../m4/common/iofile)
c
INCLUDE(common/aivb)
      character*10 charwall
      REAL fcrit
      INTEGER idetval
      REAL prlev
      common/leading/prlev
c
      ind(i,j) = max(i,j)*(max(i,j)-1)/2 + min(i,j)
c
_IF(parallel)
      logical oroot
_ENDIF
c
c...  find out the highest index of the occupied orbitals (imax)
c
_IF(atmol)
      rewind 25
      read(25) nelec,nalfa,nconf,ndets,nstruc,maxdet,maxstr,
     &         nwpack,ncoeff,imax
_ELSE
      call wr15vb(nelec,nalfa,nconf,ndets,nstruc,maxdet,maxstr,
     &            nwpack,ncoeff,imax,ndum,ndum2,'read')
_ENDIF
c
      if (nitscf.eq.maxscf) then
         write(iwr,600) nitscf,cpulft(1),charwall(),title,0.0d0
      else
         write(iwr,600) nitscf,cpulft(1),charwall(),title,evb
      end if
600   format(//,1x,27('_-_-'),
     &      //,' final vbscf results after cycle ',i5,' at ',f12.3,
     &        ' seconds ',a10,' wall',//,1x,' ** ',10a8,' **',
     &      //,' vb-energy :',f24.14,' hartree')
      if (super_cas.or.super_act) write(iwr,605)
605   format(  ' *** restrictions in the orbital optimisation ***')
c
      ki = 1
      ks = ki + nstruc
      call filnat(q(ki),nstruc)
      call submat(q(ks),nstruc,q(ki),ihfile,ihbloc,nstruc,'s')
_IF(parallel) 
_IFN(peigss)
      call pg_dgop(17193,q(ks),nstruc*(nstruc+1)/2,'+')
_ENDIF
_ENDIF
      kwgn = igmem_alloc_inf(nstruc,'vbin.m','scfina','kwgn',
     &       IGMEM_DEBUG)
      ksscr = igmem_alloc_inf(nstruc*nstruc*2,'vbin.m','scfina',
     &        'ksscr',IGMEM_DEBUG)
      ksinv = igmem_alloc_inf(nstruc*nstruc,'vbin.m','scfina',
     &        'ksinv',IGMEM_DEBUG)
      ksv1 = igmem_alloc_inf(nstruc,'vbin.m','scfina','ksv1',
     &       IGMEM_DEBUG)
      ksv2 = igmem_alloc_inf(nstruc,'vbin.m','scfina','ksv2',
     &       IGMEM_DEBUG)
      call vclr(qq(ksv2),1,nstruc)
      call submat(qq(ksscr),nstruc,q(ki),ihfile,ihbloc,nstruc,'s')
_IF(parallel) 
      call pg_dgop(17194,qq(ksscr),nstruc*(nstruc+1)/2,'+')
_ENDIF
      call square(qq(ksinv),qq(ksscr),nstruc,nstruc)

      fcrit = 10E-8
      call osinv_vb(qq(ksinv),nstruc,idetval,fcrit,qq(ksv1),qq(ksv2))

      it = 1
      psinorm = 0.0d0
      cinorm = 0.0d0
      do ij=1,nstruc
        psinorm = psinorm + ci(ij)*ci(ij)/qq(ksinv+it-1)
        cinorm = cinorm + ci(ij)*ci(ij)
        it = it + nstruc + 1
      end do

      it = 1
      do ij=1,nstruc
        qq(kwgn+ij-1) = ci(ij)*ci(ij)/qq(ksinv+it-1)
        qq(kwgn+ij-1) = qq(kwgn+ij-1)/psinorm
        it = it + nstruc + 1
      end do

      kww = ks + nstruc*(nstruc+1)/2
      do 700 i=1,nstruc
         kw = kww + i-1
         q(kw) = 0.0d0
         do 699 j=1,nstruc
            q(kw) = q(kw) + ci(i) * q(ks+ind(i,j)-1) * ci(j)
699      continue
700   continue
c
      write(iwr,608) (ci(ij),ij=1,nstruc)
608   format(/,' vb-vector :',(t13,10f10.6))
      write(iwr,771) (q(ij),ij=kww,kww+nstruc-1)
      write(iwr,772) (qq(ij),ij=kwgn,kwgn+nstruc-1)
771   format(/,' weights CC:',(t13,10f10.6))
772   format(/,' weights GN:',(t13,10f10.6))
      write(iwr,78) (q(ks+ind(i,i)-1),i=1,nstruc)
78    format(/,' struc-norm:',(t13,10f10.6))
      kbonds = kww + nstruc
      kscr = kbonds + nelec*nstruc
      kend = kscr + nelec*nstruc
      call prweig(q(kww),q(kbonds),q(kscr),nstruc,qq(kwgn))

      if ( aivb_set .eqv. .true. ) then
        write(iwr,2) prlev
2       format(/,' list contains diagrams for G-N weights larger than ',
     &        f11.8)
        write(iwr,3)
3      format(/,'  s#  |weight Ch-C |weight G-N | Atomic states ',
     & /,'------+------------+-----------+--',
     &   '------------------------------------------',
     &   '-------------------------------')
        do ijk=1,nstruc
          if ( nstruc .lt. 20 ) then
            write(iwr,'(i5,a,f10.6,a,f10.6,a,a)')
     &      ijk,'  ',q(kww+ijk-1),'  ',qq(kwgn+ijk-1),'    ',
     &      aivb_confdescr(ijk,3)
          else
            if ( qq(kwgn+ijk-1).gt.prlev) then
              write(iwr,'(i5,a,f10.6,a,f10.6,a,a)')
     &        ijk,'  ',q(kww+ijk-1),'  ',qq(kwgn+ijk-1),'    ',
     &        aivb_confdescr(ijk,3)
            end if
          end if
        end do
      end if
      call gmem_free_inf(ksv2,'vbin.m','scfina','ksv2')
      call gmem_free_inf(ksv1,'vbin.m','scfina','ksv1')
      call gmem_free_inf(ksinv,'vbin.m','scfina','ksinv')
      call gmem_free_inf(ksscr,'vbin.m','scfina','ksscr')
      call gmem_free_inf(kwgn,'vbin.m','scfina','kwgn')

      write(iwr,612) isecv
612   format(/,1x,27('-_-_'),/,1x,' vectors dumped to section ',i4)
      if (ncore.eq.0) then
         write(iwr,602)
602      format(//,' final orbitals :')
         call prvc(v,ncore+nsa,nbasis,v,'o','l')
         if ( dumpvbmo ) then
           write(iwr,'(/,a)') 
     &'   ==> dumping vb orbitals to separate file'
           write(iwr,'(a,a)') '   filename: ',vbmofile
_IF(parallel)
           if ( oroot() ) then
             call dumporb(v,ncore+nsa,nbasis,nbasis)
             goto 6021
           else
             goto 6021
           end if
_ENDIF
           call dumporb(v,ncore+nsa,nbasis,nbasis)
6021       continue
           close(vbmounit)
         end if
         write(iwr,603)
603      format(/,' orbital overlap matrix :')
      else
         write(iwr,609)
609      format(//,' final orbitals (including core) :')
         call prvc(v,ncore+nsa,nbasis,v,'o','l')
         if ( dumpvbmo ) then
           write(iwr,'(/,a)') 
     &'   ==> dumping vb orbitals to separate file (including core)'
           write(iwr,'(a,a)') '   filename: ',vbmofile
_IF(parallel)
           if ( oroot() ) then
             call dumporb(v,ncore+nsa,nbasis,nbasis)
             goto 6091
           else
             goto 6091
           end if
_ENDIF
           call dumporb(v,ncore+nsa,nbasis,nbasis)
6091       continue
           close(vbmounit)
         end if
         write(iwr,610)
610      format(/,' orbital overlap matrix (without core) :')
      end if
      call tripri(s,imax)
c
      if (lchbas(v,q(kend),ncore+nsa,'off').eq.1) then
         write(iwr,611)
611      format(//,' final orbitals (original basis) :') 
         call prvc(v,ncore+nsa,nbasis,v,'o','l')
      end if
c
cbug  dettri spoils s, which is still needed in sym1
cbug  can result in wrong natural orbitals in some cases
cbug  use a copy of s instead
c
      ks = imax + imax + 1
      kt = imax * imax + imax
      kt = kt / 2
      call fmove(s,q(ks),kt)
      call dettri(q(ks),q,q(imax+1),imax,det,irank)
      write(iwr,604) imax-irank,det
604   format(/,' nullity and determinant :',i3,2x,f7.5)
c
      kcin = igmem_alloc_inf(nstruc,'vbin.m','scfina','kcin',
     &       IGMEM_DEBUG)

c     if (nstruc.le.200) nstate = max(nstruc,nstate)
c     if (nstate.ne.0) call prhsmat(nstate,q,evb,ci,qq(kcin))
      nstate = max(nstruc,nstate)
      call prhsmat(nstate,q,evb,ci,qq(kcin))

      write(iwr,607) (qq(kcin+ij-1),ij=1,nstruc)
607   format(/,' vb-vector (now normalised): ',(t30,10f10.6))

      call gmem_free_inf(kcin,'vbin.m','scfina','kcin')
      call prnatorb(q)
c
      return
      end
      subroutine prhsmat(nstate,q,evb,ci,cin)
c
      implicit REAL (a-h,o-z), integer(i-n)
c
c...  print h and s matrices (over structures) for vbscf/vbci
c...  print eigenvalues/vectors
c
      dimension q(*),ci(*),cin(*)
c
INCLUDE(common/hsinfo)
INCLUDE(common/vbcri)
INCLUDE(common/ffile)
INCLUDE(../m4/common/iofile)
c
      ind(i,j) = max(i,j)*(max(i,j)-1)/2 + min(i,j)
c
      npri = min(nstate,nstruc)
c
      ki = 1
      ks = ki + nstruc+1
      kh = ks + nstruc*(nstruc+1)/2
      kvec = kh + nstruc*(nstruc+1)/2
      ke = kvec + nstruc*nstruc
      kscr = ke + nstruc
      ksc2 = kscr + nstruc**2
      ksc3 = ksc2 + nstruc**2
      kend = ksc3 + nstruc**2
c
      call filnat(q(ki),nstruc)
      call submat(q(ks),nstruc,q(ki),ihfile,ihbloc,nstruc,'s')
      call submat(q(kh),nstruc,q(ki),ihfile,ihbloc,nstruc,'h')
_IF(parallel)
_IFN(peigss)
      call pg_dgop(17192,q(ks),nstruc*(nstruc+1)/2,'+')
      call pg_dgop(17081,q(kh),nstruc*(nstruc+1)/2,'+')
_ENDIF
_ENDIF


c....     normalise and add nuclear repulsion to h-matrix
c..       add frozen core contribution
      do 900 i=1,nstruc
900   q(kscr+i-1) = 1/dsqrt( q(ks+i*(i+1)/2-1) )
      do 840 i=1,nstruc
         do 830 j=1,i
            rnorm = q(kscr+i-1) * q(kscr+j-1)
            ip = i*(i-1)/2 + j - 1
            q(ks+ip) = q(ks+ip) * rnorm
            q(kh+ip) = q(kh+ip) * rnorm + q(ks+ip) * core
830      continue
840   continue

      cinnorm = 0.0d0
      do i=1,nstruc
        ip = i*(i-1)/2 + i - 1
        cinnorm = cinnorm + ci(i)*ci(i)*q(ks+ip)
        do j=1,i-1
          ip = i*(i-1)/2 + j - 1
          cinnorm = cinnorm + 2*ci(i)*ci(j)*q(ks+ip)
        end do
      end do

      do i=1,nstruc
        cin(i) = ci(i)/DSQRT(cinnorm)
      end do

      znorm = 0.0d0
      do i=1,nstruc
        znorm = znorm + cin(i)*cin(i)
        do j=1,i-1
          znorm = znorm + 2*cin(i)*cin(j)*q(ks+ind(i,j)-1)
        end do
      end do

      write(iwr,711)
711   format(//,' matrix representation of the hamiltonian',
     1          ' (now normalised) :')
      write(75)nstruc
      l1=nstruc*(nstruc+1)/2
      write(75)(q(kh+i-1),i=1,l1)
      call tripri(q(kh),nstruc)
      if ( nstruc .gt. 1 ) then
        write(iwr,*) ' '
        call prenres(q(kh),nstruc,evb)
      end if
      write(iwr,722)
722   format(//,' corresponding metric :')
      write(75)(q(ks+i-1),i=1,l1)
      call tripri(q(ks),nstruc)
      call jacobs(q(kh),nstruc,q(ks),q(kvec),
     &           q(ke),2,crilow,q(kscr))
      write(iwr,701) npri
701   format(//,' lowest',i3,' eigenvectors/values :')
      call prvc(q(kvec),npri,nstruc,q(ke),'v','a')
c
c...    also print on orthogonal basis

      call submat(q(kh),nstruc,q(ki),ihfile,ihbloc,nstruc,'h')
      call submat(q(ks),nstruc,q(ki),ihfile,ihbloc,nstruc,'s')
_IF(parallel)
_IFN(peigss)
      call pg_dgop(17192,q(ks),nstruc*(nstruc+1)/2,'+')
      call pg_dgop(17081,q(kh),nstruc*(nstruc+1)/2,'+')
_ENDIF
_ENDIF
c....     add nuclear repulsion and/or frozen core to h-matrix
      do i=1,nstruc
         do  j=1,i
           ip = i*(i-1)/2 + j - 1
           q(kh+ip) = q(kh+ip) + q(ks+ip) * core
         end do
      end do
c
c        generate orthogonal matrix
c
      call extr(q(kvec),q(kvec),0,nstruc)
      call lowdins(q(kvec),q(ks),q(ksc3),q(kscr),q(ksc2),q(ke),
     1             nstruc,nstruc,crilow)
c
c...     orthogonalising transformation in q(kvec)
c...     now transform h-matrix to lowdin orthogonalised basis
c
      call dagger(nstruc,nstruc,q(kvec),nstruc,q(kscr),nstruc)
      call mult1(q(kh),q(ks),q(ksc2),nstruc,nstruc,nstruc,
     1           q(kvec),q(kscr))
c
c....   now orthogonalised h-matrix in q(ks)
c
      write(iwr,712)
712   format(//,' **(lowdin)orthogonalised representation of the ',
     1          'hamiltonian **')
      call tripr3(q(ks),nstruc)
c
      call filiky(q(ki),nstruc+1)
      call jacobt(q(ks),q(ki),nstruc,q(kvec),nstruc,q(ke),
     &            0,2,crilow,q(kscr))
      write(iwr,701) npri
      call prvc(q(kvec),npri,nstruc,q(ke),'v','a')
c
      return
      end

      subroutine prheigval(nstate,q,weightgn)
c
      implicit REAL (a-h,o-z), integer(i-n)
c
c...  print h diagonal elements for structures with weight above prlev
c
      dimension q(*)
      REAL weightgn(*)
      integer nstate,prlev
c
INCLUDE(common/hsinfo)
INCLUDE(common/vbcri)
INCLUDE(common/ffile)
INCLUDE(../m4/common/iofile)
c
      npri = min(nstate,nstruc)
c
      ki = 1
      ks = ki + nstruc+1
      kh = ks + nstruc*(nstruc+1)/2
      kvec = kh + nstruc*(nstruc+1)/2
      kscr = kvec + nstruc*nstruc
c
      call filnat(q(ki),nstruc)
      call submat(q(ks),nstruc,q(ki),ihfile,ihbloc,nstruc,'s')
      call submat(q(kh),nstruc,q(ki),ihfile,ihbloc,nstruc,'h')
_IF(parallel)
_IFN(peigss)
      call pg_dgop(17192,q(ks),nstruc*(nstruc+1)/2,'+')
      call pg_dgop(17081,q(kh),nstruc*(nstruc+1)/2,'+')
_ENDIF
_ENDIF
c...       normalise and add nuclear repulsion to h-matrix
c...       add frozen core contribution
      do 900 i=1,nstruc
900   q(kscr+i-1) = 1/dsqrt( q(ks+i*(i+1)/2-1) )
      do 840 i=1,nstruc
         do 830 j=1,i
            rnorm = q(kscr+i-1) * q(kscr+j-1)
            ip = i*(i-1)/2 + j - 1
            q(ks+ip) = q(ks+ip) * rnorm
            q(kh+ip) = q(kh+ip) * rnorm + q(ks+ip) * core
830      continue
840   continue
      write(iwr,711)
711   format(//,' diagonal elements of the hamiltonian',
     1          ' (now normalised) :')
      call prdiag(q(kh),nstruc,weightgn)
c
      return
      end


      subroutine prnatorb(q)
c
      implicit none
c
      REAL q
      dimension q(*)
      integer kvnat,keval,koccu,kscr,ndum1,ndum2,ndum3,ndum4,
     &        nsave,nocc,i
INCLUDE(common/tractlt)
INCLUDE(../m4/common/atmol3)
INCLUDE(../m4/common/iofile)
INCLUDE(../m4/common/restar)
c
      kvnat = 1
      keval = kvnat + nbasis*nbasis
      koccu = keval + nbasis
      kscr  = koccu + nbasis
      nsave = nprint
      nprint = -5
      call getq(q(kvnat),q(keval),q(koccu),ndum1,ndum2,3,
     &          ndum3,ndum4,mouta,'natorb  ')
      nprint = nsave
c      print *,'========= ndum1: ',ndum1
c      print *,'========= ndum2: ',ndum2
c      print *,'========= nbasis: ',nbasis
      do nocc=0,nbasis-1
c         print *,'================keval: ',q(keval+nocc-2)
         if (q(keval+nocc).lt.0.00000001) goto 110
      enddo
110   if (iprint.ge.10) then
         write(iwr,120)
120      format(/,' natural orbitals :')
         call prvc(q(kvnat),nocc,nbasis,q(keval),'v','l')
      else
         i=1
         write(iwr,130)
         if (nocc.le.1) then
	    write(iwr,140)
	    write(iwr,145)
            write(iwr,150) 1,q(keval)
	 else if (nocc.le.2) then
	    write(iwr,160)
	    write(iwr,165)
            write(iwr,170) 1,q(keval),2,q(keval+1)
	 else
	    write(iwr,180)
	    write(iwr,185)
125         write(iwr,190) i,q(keval+i-1),i+1,
     &         	  q(keval+i),i+2,q(keval+i+1)
            i=i+3
            if ((nocc-i).ge.2) goto 125
	    if ((nocc-i).eq.1) then
               write(iwr,170) i,q(keval+i-1),i+1,q(keval+i)
	    end if
	    if ((nocc-i).eq.0) then
               write(iwr,150) i,q(keval+i-1)
	    end if
	 end if
130   format(/,' natural orbital occupancies :')
140   format(/,' orb | occupancy |')
145   format(' ----+-----------+')
150   format(i4,' |',f9.4,'  |')
160   format(/,' orb | occupancy | orb | occupancy |')
165   format(' ----+-----------+-----+-----------+')
170   format(i4,' |',f9.4,'  |',i4,' |',f9.4,'  |')
180   format(/,' orb | occupancy | orb | occupancy | orb',
     &           ' | occupancy |')
185   format(' ----+-----------+-----+-----------+----',
     &           '-+-----------+')
190   format(i4,' |',f9.4,'  |',i4,' |',f9.4,'  |',i4,' |',f9.4,'  |')
      endif
      return
      end

      subroutine dettri(s,ir,ic,ndim,det,irank)
      implicit REAL  (a-h,o-z) ,  integer   (i-n)
      dimension s(*),ir(*),ic(*)
      parameter (crit=1.0d-3)
      ind(i,j) = max(i,j) * (max(i,j)-1) / 2 + min(i,j)
      ipar = 1
      irank = 0
      do 1 i=1,ndim
         ir(i) = i
1     ic(i) = i
      det = 1.0d0
      do 50 i=1,ndim
         call pivtri(s,ir,ic,ndim,i,ipar)
         diag = s(ind(i,i))
         det = det * diag
         if (dabs(det).lt.crit) then
            det = det * ipar
            return
         end if
         irank = irank + 1
         do 30 j=i+1,ndim
            gauss = -s(ind(i,j))/diag
            s(ind(j,j)) = s(ind(j,j)) + gauss * s(ind(i,j))
            do 20 k=j+1,ndim
               s(ind(k,j)) = s(ind(k,j)) + s(ind(k,i))*gauss
20          continue
30       continue
50    continue
      det = det * ipar
      return
      end
      subroutine bond1(iscr,lword)
      implicit REAL  (a-h,o-z) ,  integer   (i-n)
INCLUDE(common/c8_16vb)
c.....
c.....find out which singly occupied orbitals are bonded via spin-
c.....functions.
c.....bond1 calls bond2 as a remedy against headache
c.....
      dimension iscr(lword)
      call getpa2 (nelec,nalfa,iscr,npack,nconf,nstruc,ndett,lword,imax)
c..   # dets. per structure
      kndps = 1
c..   det.numbers per structure
      kidps = kndps + nstruc
c..   # dets. per configuration
      kndpc = kidps + ndett
c..   # structures per configuration
      knspc = kndpc + nconf
c..   packed dets.
      kpdet = knspc + nconf
c..   unpacked det (1)
      kdet = kpdet + npack
c..   singly occ. orbitals
      ksing = kdet + nelec
      n1 = 0
      ip = kpdet
      do 20 i=kndpc,kndpc+nconf-1
         call izero(nelec,iscr(kdet),1)
         call unpack(iscr(ip),n8_16,iscr(kdet),nelec)
capollo  mind 64 bits per packed word => nipw()
         ip = ip + (((iscr(i)*nelec-1)/(64/n8_16))+1)*nipw()
         do 10 j=1,imax
            nj = ioccup(j,iscr(kdet),nelec)
            if (nj.eq.1) then
               ibefor = locati(j,iscr(ksing),n1)
               if (ibefor.eq.0) then
                  iscr(ksing+n1) = j
                  n1 = n1 + 1
               end if
            end if
10       continue
20    continue
      call icopy(n1,iscr(ksing),1,iscr(kdet),1)
      ksing = kdet
      ksc = ksing + n1
      llwor = lword-ksc
      call bond2(nconf,nelec,nalfa,iscr(kndps),iscr(kndpc),iscr(knspc),
     &          iscr(kidps),iscr(kpdet),iscr(ksing),n1,iscr(ksc),llwor)
      return
      end
      subroutine bond2(nconf,nelec,nalfa,ndetps,ndetpc,nstrpc,idetps,
     &                 pdet,ione,none,idet,lword)
      implicit REAL  (a-h,o-z) ,  integer   (i-n)
      dimension ndetps(*),ndetpc(*),nstrpc(*),pdet(*),ione(*),
     &          idet(nelec,*),idetps(*)
INCLUDE(../m4/common/sizes)
INCLUDE(../m4/common/iofile)
INCLUDE(common/turtleparam)
INCLUDE(common/c8_16vb)
INCLUDE(common/vbbonds)
INCLUDE(common/twice)
INCLUDE(common/infato)
      logical samesp
      samesp(i,j) = ((i.le.nalfa.and.j.le.nalfa).or.
     &               (i.gt.nalfa.and.j.gt.nalfa)    )
      if (hguess) then
c
c..      fill isingly for subroutine hybmak
c..      this is undone by hybmak
c
         call icopy(none,ione,1,isingly,1)
         nsingly = none
      end if
c.....
      if (nbonds.ne.0) go to 80
      nbonds = 0
      do 70 i=1,none-1
c        loop over singly occs.
         ifirst = nbonds+1
         nbondi = 0
         i1 = ione(i)
         do 60 j=i+1,none
            i2 = ione(j)
c           loop over configs.
            it = 0
            ip = 1
            im = 0
            do 50 k=1,nconf
               if (nelec*ndetpc(k).gt.lword)
     &                    call corfait(nelec*ndetpc(k),lword,'bond2')
               call izero(nelec*ndetpc(k),idet,1)
               call unpack(pdet(ip),n8_16,idet,nelec*ndetpc(k))
               ip = ip + ( (ndetpc(k)*nelec-1)/(64/n8_16) ) + 1
c              loop over structures per config.
               do 40 l=1,nstrpc(k)
c                 loop over dets. per structure
                  it = it + 1
                  do 30 mm=1,ndetps(it)
                     m = idetps(im+mm)
                     nsing1 = ioccup(i1,idet(1,m),nelec)
                     nsing2 = ioccup(i2,idet(1,m),nelec)
                     ipos1  = locati(i1,idet(1,m),nelec)
                     ipos2  = locati(i2,idet(1,m),nelec)
                     if (samesp(ipos1,ipos2).or.
     &                   nsing1.ne.1.or.nsing2.ne.1) then
c....                   i1 and i2 are not bonded in this structure
                        im = im + ndetps(it)
                        goto 40
                     end if
30                continue
                  im = im + ndetps(it)
                  jbefor = locati(j,ibond(ifirst,2),nbondi)
                  if (jbefor.eq.0) then
                     nbonds = nbonds + 1
                     if (nbonds.gt.mbonds)call vberr(' too many bonds')
                     nbondi = nbondi + 1
                     ibond(nbonds,1) = i1
                     ibond(nbonds,2) = i2
                     goto 60
                  end if
40             continue
50          continue
60       continue
70    continue
c...  print results
80    if (nbonds.gt.0) 
     1    write(iwr,22) (ibond(ij,1),'-',ibond(ij,2),ij=1,nbonds)
22    format(/,' spin-bonds :',(t14,15(1x,i2,a1,i2,1x)),/)
      return
      end

c***********************************************************************
      subroutine clvecvb(v,nv,nbasis,set,test,caller)
c
c...  clear the part of the vector that does not belong to th hybrid
c...  also for active vectors; check if hybrid definition is OK
c...  if test='small' only minor deviations are acceptable
c
      implicit REAL (a-h,o-z), integer (i-n)
INCLUDE(../m4/common/iofile)
INCLUDE(../m4/common/sizes)
INCLUDE(common/turtleparam)
INCLUDE(common/infato)
INCLUDE(common/twice)
INCLUDE(common/vbcri)
      common /n2hybr/ n2hyb
      dimension v(nbasis,nv)
      dimension iat(4),jat(4)
      character*(*) set,test,caller
      character*80 text
      logical ochatao,oflop
c
      if (n2hyb.gt.1) then
        write(iwr,'(a,I3,a)') ' WARNING!! AO occurs in ',n2hyb,
     &                        ' hybrid definitions.'
        write(iwr,*) 
      end if
      if (test.eq.'small') then
         small = crihyb*100.0d0
      else if (test.eq.'null') then
         small = 0.0d0
      else if (test.eq.'any') then
         small = 1.0d99
      else
         call caserr('wrong call to clvecvb')
      end if
      if (.not.hybry) call caserr('chhybrid called for non-hybry')
      toobig = small*1000.0d0
      amm = -1.0d0
c    
      if (.not.hybry) call caserr('clvecvb call for non-hybrid')
c.....
c.....clear delocalised orbitals so that they have non-zero coeffs.
c.....on their atoms only or check if hybrids are ok
c.....
      do 30 i=1,nv
c.....   find largest coefficient
         call abmax(v(1,i),nbasis,1,'sq',1,am,im)
c...       if active set only is done
         if (set.eq.'active') then
cremco           oflop = ochatao(im,i,'iacat')
            if (ofrzmo(i)) go to 30
         end if
c
         do j=1,nbasis
            if (.not.ochatao(im,j,'common')) then
               if (dabs(v(j,i)).gt.amm) then
                  amm = dabs(v(j,i))
                  jmm = j
                  imm = i
                  immm = im
               end if
               if (set.ne.'check') v(j,i) = 0.0d0
            end if
         end do
c...     see if we accidently deleted the complete vector
            dn = ddot(nbasis,v(1,i),1,v(1,i),1)
            if (dn.lt.1.0d-10) then
               write(iwr,*) ' *** in cleaning mo ',i,
     1                      ' the mo is gone - norm ',dn,' ***'
               write(iwr,*) ' *** CHECK orbital/hybrid definitions ***'
               call vberr(' clvecvb deleted orbital ')
            end if
c
30     continue
c
      if (amm.gt.small) then
         write(text,'(a14,i4,a14,i4,a4,1pd12.5,a8,i4,2x,a12)')   
     1         ' hybrid-vector',imm,' with major ao',immm,
     2         ' has',amm,' at pos.',jmm,caller
         write(iwr,'(1x,80(1h*),/,a80)') text
            oflop = ochatao(immm,immm,'print')
            oflop = ochatao(jmm,jmm,'print')
         write(iwr,'(1x,80(1h*))') 
         if (amm.gt.toobig) then
            call prvc(v,nv,nbasis,v,'o','l')
            call caserr(text)
         end if
      end if
c
      return
      end
c***********************************************************************
      logical function ochatao(iin,jin,op)
c
c     op = 'set'   : set atomao array
c     op = 'iacat' : check if mo jin is in the iacat of atoms for iin
c     op = 'atomao': check if (one of the) atom(s) of jin is iin 
c...  op = 'common': check of ao's iin and jin are in common set 
c...  op = 'print' : print the hybrid affilations of ao's iin and jin 
c...  used in clvecvb,clmetr : required because hybrids may now overlap
c...           (atomao(j).eq.atomao(j)) 
c
      implicit REAL (a-h,p-z), integer (i-n), logical (o)
c
      integer iin,jin
      character*(*) op
c
INCLUDE(../m4/common/sizes)
INCLUDE(common/turtleparam)
INCLUDE(common/infato)
INCLUDE(common/tractlt)
INCLUDE(../m4/common/iofile)
      common /n2hybr/ n2hyb
      dimension iat(9),jat(9)
c
      if (op.eq.'set') then
         call izero(nbasis,atomao,1)
         n2hyb = 0
         do i=1,natom
            nnn = 0
            do j=1,nopa(i)
               if (atomao(iopa(j,i)).ne.0) then
                  nn = atomao(iopa(j,i))
                  nnn = 0
1                 nn = nn / maxato
                  nnn = nnn + 1
                  if (nn.gt.0) go to 1
                  atomao(iopa(j,i)) = atomao(iopa(j,i)) + i*maxato**nnn
               else
                  atomao(iopa(j,i)) = i
               end if
               n2hyb = max0(n2hyb,nnn)
            end do
         end do
         n2hyb = n2hyb + 1
         if (n2hyb.gt.1) then
            write(iwr,6905) n2hyb
6905        format(1x,/,1x,45(1h*),/,
     1             ' * Some aos occur in',i6,' hybrids           *',/,
     2             ' * This can yield inconsistent results       *',/,
     2             ' * The SUPER option is recommended           *',/,
     4             1x,45(1h*),/,1x)
         end if
      else if (op.eq.'print') then
         do k=iin,jin
            iatom = atomao(k)
            nni = 0
5           nni = nni + 1
            if (nni.gt.9) call caserr(' overflow in ochatao ')
            iat(nni) = mod(iatom,maxato)
            iatom = iatom/maxato
            if (iatom.gt.0) go to 5
            write(iwr,600) k,(iat(kk),atoms(iat(kk)),kk=1,nni)
600         format(' AO ',i5,' hybrids ',9(i5,1x,a8))
         end do
      else if (op.eq.'iacat') then
         iatom = atomao(iin)
         nni = 0
10       nni = nni + 1
            if (nni.gt.9) call caserr(' overflow in ochatao ')
            iat(nni) = mod(iatom,maxato)
            iatom = iatom/maxato
         if (iatom.gt.0) go to 10
         ithere = 0
         do k=1,nni
            ithere=ithere+locati(jin,iacat(1,iat(k)),nacat(iat(k)))
         end do
         if (ithere.eq.0) call caserr('hybrid foulup')
      else if (op.eq.'atomao'.or.op.eq.'common') then
         if (op.eq.'atomao') then
            iat(1) = iin
            nni = 1
         else        
            iatom = atomao(iin)
            nni = 0
20          nni = nni + 1
               if (nni.gt.9) call caserr(' overflow in ochatao ')
               iat(nni) = mod(iatom,maxato)
               iatom = iatom/maxato
            if (iatom.gt.0) go to 20
         end if
         jatom = atomao(jin)
         nnj = 0
30       nnj = nnj + 1
            if (nnj.gt.9) call caserr(' overflow in ochatao ')
            jat(nnj) = mod(jatom,maxato)
            jatom = jatom/maxato
         if (jatom.gt.0) go to 30
c
         oo = .false.
         do ii = 1,nni
            do jj = 1,nnj
               oo = oo.or.(iat(ii).eq.jat(jj))
            end do
         end do
         ochatao = oo
      end if
c
      return
      end
      subroutine clmetr(sao,nbasis,test)
      implicit REAL (a-h,o-z), integer(i-n)
c.....
c.....strips ao-metric (remove overlaps between atoms)
c.... *  if clean is not set , maybe ignore the small tests ??
c.....
INCLUDE(../m4/common/sizes) 
INCLUDE(common/turtleparam)
INCLUDE(common/infato)
INCLUDE(common/vbcri)
      dimension sao(*)
      character*(*) test
      character*80 text
      logical ochatao
c    
      if (.not.hybry) call caserr('clmetr for non hybry')
c
      kk = 0
      do i=1,nbasis
         do j=1,i
            kk = kk + 1
            if (.not.ochatao(i,j,'common')) then
               if (test.eq.'small'.and.dabs(sao(kk)).gt.crihyb) then
                  write(text,*) 'remco: hybrids',i,j,' have',sao(kk)
                  call caserr(text)
               end if
               sao(kk) = 0.0d0
            end if
         end do
      end do
c
      return
      end
      
      subroutine doubly(iscr,lword)
      implicit REAL  (a-h,o-z),  integer  (i-n)
c.....
c.....this routine finds out about the restrictions to the use of the
c.....total function space for the improvement of the occupied orbitals,
c.....that are dictated by the pauli-principle. (together with "doubl2",
c.....see that one too)
c.....inforb contains per scf-orbital the orbitals that should be
c.....removed from the virtual space. in practice (see "redorb") the
c.....virtuals having the largest overlap with the mentioned orbitals
c.....are excluded.
c.....
INCLUDE(../m4/common/sizes)
INCLUDE(common/turtleparam)
INCLUDE(common/c8_16vb)
INCLUDE(common/scftvb)
INCLUDE(common/brill)
INCLUDE(common/twice)
      logical there
      dimension iscr(lword)
c
c..   ioror contains per orbital the orbitals it occurs with in
c..   the configurations (is used in redorb, if .super.)
c
      call izero(maxact*maxact,ioror,1)
      call getpac(nelec,nalfa,iscr,npack,ngroup,lword,imaxor,maxdet)
      if (imaxor.gt.5*maxact) call vberr('occup. orbital index > 100')
      ip   = ngroup + 1
      idet = ip + npack
      nlast = idet + nelec*maxdet
      if (nlast.gt.lword) call corfait(nlast,lword,'in doubly')
      ione = 1
      call setsto(maxact*maxact*5,ione,inforb)
      ndoubly = imaxor
      nsingly = 0
      do 1 i=1,ndoubly
1     idoubly(i) = i
c.....idoubly will contain the doubly occupieds (see sets)
c.....
c.....partitioning iscr
c..... dets per configuration     iscr(1)                    ngroup long
c..... packed dets.               iscr(ngroup+1)             npack   ,,
c..... unpacked det               iscr(ngroup+npack+1)       nelec   ,,
c.....
      do 40 i=1,ngroup
c.....
         call izero(nelec*iscr(i),iscr(idet),1)
         call unpack(iscr(ip),n8_16,iscr(idet),nelec*iscr(i))
c
capollo mind 64 bits per packed word => nipw()
         ip = ip + ((iscr(i)*nelec-1) / (64/n8_16) + 1) * nipw()
c..      fill ioror
         do 43 k=1,nelec-1
            korb = iscr(idet+k-1)
            do 42 l=k+1,nelec
               lorb = iscr(idet+l-1)
               ioror(korb,lorb) = 1
               ioror(lorb,korb) = 1
42          continue
43       continue
         do 30 k=1,iscr(i)
            kk = (k-1)*nelec
            do 20 j=1,imaxor
               ioccur = ioccup(j,iscr(idet+kk),nelec)
               if ( ioccur.lt.2) then
c.....            store this scf-orbital as having single or variable
c.....            occupancy. this means it will be used as a virtual
c.....            orbital
                  if (locati(j,isingly,nsingly).eq.0) then
                     nsingly = nsingly + 1
                     isingly(nsingly) = j
                  end if
c......           if it was a doubly occupied orbital delete it from
c.....            this list
                  idel = locati(j,idoubly,ndoubly)
                   if (idel.ne.0) then
                     ndel = ndoubly - idel
                     ndoubly = ndoubly - 1
                     call icopy(ndel,idoubly(idel+1),1,idoubly(idel),1)
                  end if
               end if
               if (ioccur.le.1) then
c.....            undo restrictions involving this orbital for all
c.....            other orbitals in this configuration
                  do 10 l=idet,idet+nelec-1
                     if (iscr(l+kk).ne.j) then
c.....                  find the position of the scf orbital
                        inforb(j,iscr(l+kk)) = 0
                     end if
10                continue
               end if
20          continue
30       continue
40    continue
      return
c.....this job will be completed by "doubl2"
      end
      subroutine doubl2
      implicit REAL  (a-h,o-z),  integer  (i-n)
c.....
c..... finishes the job of doubly
c.....
INCLUDE(../m4/common/sizes)
INCLUDE(../m4/common/iofile)
INCLUDE(common/turtleparam)
INCLUDE(common/c8_16vb)
INCLUDE(common/scftvb)
INCLUDE(common/brill)
INCLUDE(common/twice)
      logical there
c.....
c.....condense information
c.....
      i = 0
      do 100 ii=1,nequi
         istart = i
         do 50 jj=1,iequi(ii)
            i = i + 1
            ninfor(i) = 0
            do 40 j=1,nscf
               if (inforb(j,i).ne.0) then
                  if (iexw(1,i).ne.j) then
                     ninfor(i) = ninfor(i) + 1
                     inforb(ninfor(i),i) = j
                  end if
               end if
40          continue
            if (reduce) then
               ninfor(i) = ninfor(i) + 1
               inforb(ninfor(i),i) = iexw(1,i)
            end if
50       continue
         minres = ninfor(istart+1)
         minpos = istart + 1
         do 60 k=istart+2,istart+iequi(ii)
            if (ninfor(k).lt.minres) then
               minpos = k
               minres = ninfor(k)
            end if
60       continue
         if (iequi(ii).gt.1) then
c.....      the freeest orbital liberates the others
c.....      locate freeest orbital on inteq and translate the
c.....      restrictions (simple horizontal correspondence)
            ifree = 0
            do 65 l=1,ninter
65          if (inteq(1,l).eq.iexw(1,minpos)) ifree = l
            if (ifree.eq.0) then
               write(iwr,44) iexw(1,minpos)
44             format(' cant find orbital',i3,' in internal',
     &                ' equivalence definitions')
               call vberr('please complete equivalence information')
            end if
            do 90 k=1,ninfor(minpos)
c....          locate the k-th restriction on inteq(*,ifree)
               ipos1 = 0
               do 70 l=2,ninteq(ifree)+1
70             if (inteq(l,ifree).eq.inforb(k,minpos)) ipos1 = l
               if (ipos1.eq.0) then
c....             forget about it, this excitation should not occur,
c....             see "equiv"
                  goto 100
               end if
               do 80 l=istart+1,istart+iequi(ii)
                  if (l.ne.minpos) then
c....                locate l-th orbital on inteq and use ipos2
                     ipos2 = 0
                     do 75 m=1,ninter
75                   if (inteq(1,m).eq.iexw(1,l)) ipos2 = m
                     if (ipos2.eq.0) then
c....
c....                   forget about it, this excitation should not
c....                   occur, see "equiv"
c....
                        goto 100
                     end if
                     inforb(k,l) = inteq(ipos1,ipos2)
                     ninfor(l) = ninfor(minpos)
                  end if
80             continue
90          continue
         end if
100   continue
      it = 1
      do 200 i=1,nequi
        do 190 j=1,iequi(i)
1234       nn = ninfor(it+j-1)
           do 180 k=1,nn
              if (inforb(k,it+j-1).eq.0) then
                 do 170 l=it,it+iequi(i)-1
                    ninfor(l) = ninfor(l)-1
170              call icopy(nn-k,inforb(k+1,l),1,inforb(k,l),1)
                 goto 1234
              end if
180        continue
190     continue
        it = it + iequi(i)
200   continue
      ipar = 1
      call ibubbl(idoubly,ndoubly,ipar)
      call ibubbl(isingly,nsingly,ipar)
      if (iprins.gt.1000) then
         write(iwr,11)
11       format(/' information on pauli-principle related restrictions',
     &           ' to excitations  ')
         write(iwr,22)
22       format(/' scf-orb. ,  orbitals to be removed =>')
         do 110 i=1,nscf
            if(ninfor(i).ne.0) then
               write(iwr,33) iexw(1,i),(inforb(j,i),j=1,ninfor(i))
33             format(5x,i3,5x,35i3)
            end if
110      continue
         if (ndoubly.gt.0) then
             write(iwr,55) (idoubly(ij),ij=1,ndoubly)
55           format(1x,'orbitals that are occupied twice always  ',
     &              (t45,29i3))
         end if
         if (nsingly.gt.0) then
             write(iwr,66) (isingly(ij),ij=1,nsingly)
66           format(1x,'            the other occupied orbitals  ',
     &              (t45,29i3))
         end if
         write(iwr,'(a)')' occurance of orbitals per orbital:'
         do 210 i=1,nscf
            write(iwr,77) (ioror(i,j),j=1,nscf)
77          format(1x,(35i3))
210      continue
      end if
      return
      end
      subroutine elimin(n,map,ifile,ibloch,nstruc)
      implicit REAL  (a-h,o-z) ,  integer   (i-n)
      dimension map(n)
c
c.....
c.....it is decided the csf's' indicated on map, are not such a good
c.....idea after all (coupling between variational freedom with respect
c.....to mo-coefficients and ci-coefficients in the mcscf wavefunction)
c.....delete them effectively by making the corresponding h-matrix
c.....and s-matrix elements zero. the s-diagonal elements will be put to
c.....1. (to prevent division problems, it is irrelevant otherwise)
c
c..... version for labeled matrix elements
c..... must be made more efficient (FOKKE)
c.....
c
      common/matrixt/hof(170),sof(170),icof(170),irof(170),no,iblockk
      logical och
_IF(linux)
      external fget
_ENDIF 
c
      iblock = ibloch
      call search(iblock,ifile)
10    call fget(hof,nw,ifile)
      och = .false.
      if (nw.eq.0) then
        return
      end if
      do 20 i=1,no
         do 15 j=1,n
            if (map(j).eq.irof(i).or.map(j).eq.icof(i)) then
               hof(i) = 0.0d0
               sof(i) = 0.0d0
               if (icof(i).eq.irof(i)) sof(i) = 1.0d0
               och = .true.
               go to 20
            end if
15       continue
20    continue
      if (och) then
         call search(iblock,ifile)
         call put(hof,511,ifile)
      end if
      iblock = iblock + 1
      go to 10
c
      return
      end
      subroutine equivt(vec,haomo,nao,nmo,qq)
c
c.....define internal equivalences and check phases
c
      implicit REAL (a-h,o-z), integer(i-n)
      common /fast/ nospin,nosymm,nosymc
      logical nospin,nosymm,nosymc
INCLUDE(../m4/common/sizes)
INCLUDE(../m4/common/iofile)
INCLUDE(../m4/common/gmempara)
INCLUDE(common/turtleparam)
INCLUDE(common/brill)
INCLUDE(common/scftvb)
INCLUDE(common/twice)
INCLUDE(common/vbcri)
      dimension vec(nao,*),haomo(*),qq(*)
      ind(i,j)=max(i,j)*(max(i,j)-1)/2+min(i,j)
c
c     call doubly(scr,lword) -- see comment in vbin --
c
      call toetoe(vec,haomo,nao,nmo,qq)
c
c      call tripri(haomo,nao+nmo)
c.....if toetoe is adapted so that it automaticly reorders vectors,
c.....then haomo needs to be recalculated now (!)
c
c.....make a copy of part of haomo
      n=2*nmo
c     kscr  = igmem_alloc_inf((n*n+n)/2,'vbin.m','equivt',
c    @                         'kscr',IGMEM_DEBUG)
      kscr=1
      call fmove(haomo,qq(kscr),(n*n+n)/2)
      if (super) then
         do 2 i=1,nmo
            do 1 j=1,nmo
               qq(kscr-1+ind(i+nmo,j))=qq(kscr-1+ind(i,j))
1           continue
2        continue
      end if
c
c.....define internal equivalences (inteq)
c
      do 10 i=1,nmo
         inteq(1,i)=i
10    continue

      iorb1=1
      do 50 i=1,nequi
         do 40 j=1,nmo
            lstart=j+nmo+1
            h1=dabs(qq(kscr-1+ind(j+nmo,iorb1)))
            if (h1.le.remcri) then
c.....         skip excitation iorb1 => j
               goto 40
            end if
            inteq(j+1,iorb1)=j
            do 30 k=2,iequi(i)
               iorb2=iorb1+k-1
               do 20 l=lstart,lstart+nmo-1
                  if (l.le.2*nmo) then
                     m=l
                  else
                     m=l-nmo
                  end if
                  h2=dabs(qq(kscr-1+ind(m,iorb2)))
                  if (dabs(1.0-h2/h1).lt.symcri) then
c.....               found match
                     lstart=m+1
                     inteq(j+1,iorb2)=m-nmo
                     qq(kscr-1+ind(m,iorb2))=0.0d0
                     goto 30
                  end if
20             continue
30          continue
40       continue
         iorb1=iorb1+iequi(i)
50    continue
c
c     call gmem_free_inf(kscr,'vbin.m','equivt','kscr')
c
      if (iprins.gt.9999) then
         write(iwr,60)
60       format(/,' internal equivalence definitions',/,
     &                           2x,'mo',3x,'all mos =>')
         do 70 j=1,nscf
70       write(iwr,80) (inteq(i,j),i=1,nscf+1)
80       format(1x,i3,2x,100i3)
      end if
c
c.....check phases
c
      do 90 i=1,nscf
         ieqsig(1,i)=1
90    continue
      i=1
      nam=nao+nmo
      do 140 j=1,nequi
         do 130 k=2,iequi(j)
            do 120 l=i+2,nscf+1
               if (inteq(1,i).ne.0) then
                  iorb1a=inteq(1,i)
                  iorb2a=inteq(1,i+k-1)
                  iorb1b=inteq(l,i)
                  iorb2b=inteq(l,i+k-1)
                  h1=haomo(ind(iorb1a,iorb1b))
                  h2=haomo(ind(iorb2a,iorb2b))
                  if (h1*h2.lt.0.0d0.and.iorb2b.ne.0) then
                     if (ieqsig(1,iorb2b).ne.1) then
c.....                  phase was already altered for another one
                        call vberr('inconsistent equivalence')
                     else
c.....                  change phase of iorb2b
                        ieqsig(1,iorb2b)=-ieqsig(1,iorb2b)
                        do 100 n=1,nao
                           vec(n,iorb2b)=-vec(n,iorb2b)
100                     continue
                        do 110 n=1,nam
                           if (n.ne.iorb2b)
     &                     haomo(ind(iorb2b,n))=-haomo(ind(iorb2b,n))
110                     continue
                     end if
                  end if
               end if
120         continue
130      continue
         i=i+iequi(j)
140   continue
      i=0
      do 170 j=1,nequi
         do 160 k=1,iequi(j)
            i=i+1
            do 150 l=2,nexw(i)+1
               h1=haomo(ind(iexw(1,i),iexw(l,i)))
               ieqsig(l,i)=-int(dsign(1.0d0,h1))
150         continue
160      continue
170   continue
      ninter=nscf
      do 180 i=1,nscf
         ninteq(i)=nscf
180   continue
c
      if (nosymm) then
         do 300 i=1,nscf
            exact(i)=.false.
300      continue
      end if
c
c.....finish the job of doubly
      call doubl2
c
      return
      end
      subroutine extr(vco,vcn,nold,new)
c
      implicit REAL (a-h,o-z), integer (i-n)
c
c...  add unit matrix to set of old vectors / vcn may be vco
c
      dimension vco(nold,nold),vcn(new,new)
c
      if (new.le.nold) return
      nold1 = nold + 1
c...  copy vco to vcn
      do 20 i=1,nold
         do 10 j=1,nold
            vcn(nold1-j,nold1-i) = vco(nold1-j,nold1-i)
10       continue
20    continue
c...  add unit matrix
      do 40 i=nold1,new
         do 30 j=1,new
            vcn(i,j) = 0.0d0
            vcn(j,i) = 0.0d0
30       continue
40    vcn(i,i) = 1.0d0
c
      return
      end
      subroutine filnat(if,n)
      implicit REAL (a-h,o-z),integer(i-n)
      dimension if(*)
      do 10 i=1,n
10    if(i) = i
      return
      end
      subroutine filiky(iky,n)
      implicit REAL (a-h,o-z), integer (i-n)
c.....fill iky, for triangles
      dimension iky(n)
      do 10 i=1,n
10    iky(i) = i*(i-1)/2
      return
      end
      subroutine hybmak(v,v2,ir,ic,nbasis,nsa,qq)
      implicit REAL (a-h,o-z), integer(i-n)
INCLUDE(../m4/common/sizes)
INCLUDE(../m4/common/iofile)
INCLUDE(common/turtleparam)
INCLUDE(common/brill)
INCLUDE(common/vbbonds)
      integer ibond_save(mbonds,2)
INCLUDE(common/infato)
INCLUDE(common/twice)
INCLUDE(common/scftvb)
INCLUDE(common/vbcri)
INCLUDE(../m4/common/gmemdata)
INCLUDE(../m4/common/gmempara)
c
      dimension v(nbasis,*),v2(nbasis,*),ir(nbasis),ic(nbasis),qq(*)
      itri(i,j) = max(i,j)*(max(i,j)-1)/2 + min(j,i)
c
      ksao = igmem_alloc_inf(nbasis*(nbasis+1)/2,'vbin.m','hybmak',
     &                       'ksao',IGMEM_DEBUG)
      call get1e(qq(ksao),dummy,'s',dummy)
      kvt = igmem_alloc_inf(nbasis,'vbin.m','hybmak',
     &                         'kvt',IGMEM_DEBUG)
      do i=1,nbonds
         ibond_save(i,1) = ibond(i,1)
         ibond_save(i,2) = ibond(i,2)
      end do
      if (ichhyb.gt.0) then
c....    atoms were input / find the most likely binding mo'st
c....    ir is abused for admin
         do i=1,nbasis
            ir(i) = 0
         end do
         do i=1,nbonds
            call vclr(v2,1,nbasis*nbasis)
            iatom1 = ibond(i,1)
            iatom2 = ibond(i,2)
            nat1 = nacat(iatom1)
            nat2 = nacat(iatom2)
c...  do not copy the ones already done
            do k=1,nat1
               if (ir(iacat(k,iatom1)).eq.0)
     1            call fmove(v(1,iacat(k,iatom1)),v2(1,k),nbasis)
            end do
            do k=1,nat2
               if (ir(iacat(k,iatom2)).eq.0)
     1            call fmove(v(1,iacat(k,iatom2)),v2(1,k+nat1),nbasis)
            end do
c
            ksmo = igmem_alloc_inf(nbasis*(nbasis+1)/2,'vbin.m',
     &                             'hybmak','ksmo',IGMEM_DEBUG)
            call fmos(qq(ksmo),qq(ksao),v2,qq(kvt),nat1+nat2,nbasis,
     1                 crilow)
c...  look for biggest overlap between iatom1 and iaetom2
            amax = -1.0d0
            do k=1,nat1
               do l=1,nat2
                  big = abs(qq(ksmo-1+itri(k,l+nat1)))
                  if (big.gt.amax) then
                     amax = big
                     kbig = k
                     lbig = l
                  end if
               end do
            end do
            call gmem_free_inf(ksmo,'vbin.m','hybmak','ksmo')
            ibond(i,1) = iacat(kbig,iatom1)
            ibond(i,2) = iacat(lbig,iatom2)
            ir(ibond(i,1)) = 1
            ir(ibond(i,2)) = 1
         end do
c... we emulated bonding pairs
      end if
c.....v2 is possibly a copy of v
      if (guesao) then
c...  linear combinations of raw ao's
         call vclr(v2,1,nbasis*nbasis)
         do i=1,nbasis
            v2(i,i) = 1.0d0
         end do
      else
c...  linear combinations of startvectors will be used for the
c...  guessing of the hybrids
        call fmove(v,v2,nbasis*nsa)
      end if
      do 100 i=1,nbonds
c.....
c.....   find the atoms the bond connects
c.....
            iatom1 = 0
            iatom2 = 0
            ib1 = ibond(i,1)
            ib2 = ibond(i,2)
            ibefo1 = locati(ib1,ibond(1,1),i-1)
            ibefo2 = locati(ib2,ibond(1,2),i-1)
c         if (ibefo1.ne.0.or.ibefo2.ne.0) then
c            write(iwr,55) ib1,ib2
c55          format(1x,i3,'--',i3,' is not used in hybrid construction')
c....       these hybrids are not made, so temporarily exclude them
c....       from the search-list
c            ibond(i,1) = -ibond(i,1)
c            ibond(i,2) = -ibond(i,2)
c            goto 100
c         end if
            do 20 j=1,natom
               do 10 k=1,nacat(j)
                  if (iacat(k,j).eq.ib1) iatom1 = j
                  if (iacat(k,j).eq.ib2) iatom2 = j
10             continue
20          continue
            if (iatom1.eq.0.or.iatom2.eq.0) then
               call vberr('cannot find active bonded orbitals')
            end if
c....
c....    now construct (non-hermitian) partial metric
c....
         nr = 0
         nc = 0
         if (.not.guesao) then
c....       only singly occupied orbitals are allowed to mix
            do 30 j=1,nacat(iatom1)
               iorb = iacat(j,iatom1)
               ithere = locati(iorb,isingly,nsingly)
               if (ithere.ne.0) then
                  nr = nr + 1
                  ir(nr) = iorb
                  if (iorb.eq.ib1) iorig = nr
               end if
30          continue
            do 40 j=1,nacat(iatom2)
               jorb = iacat(j,iatom2)
               ithere = locati(jorb,isingly,nsingly)
               if (ithere.ne.0) then
                  nc = nc + 1
                  ic(nc) = jorb
                  if (jorb.eq.ib2) jorig = nc
               end if
40          continue
            if (jorig.eq.0.or.iorig.eq.0) call vberr('hybmak ?')
         else
c....
c....       all ao's of the atom are used
c....  nr = nopa(iatom1)
            nc = nopa(iatom2)
            call icopy(nr,iopa(1,iatom1),1,ir,1)
            call icopy(nc,iopa(1,iatom2),1,ic,1)
            call abmax(v(1,ib1),nbasis,1,'sq',1,amaxx,imaxx)
            iorig = imaxx
            call abmax(v(1,ib2),nbasis,1,'sq',1,amaxx,imaxx)
            jorig = imaxx
         end if
c
         ksvd = igmem_alloc_inf(nc*nr,'vbin.m','hybmak',
     &                          'ksvd',IGMEM_DEBUG)
c....
c....    always let nc be the least
c....
         if (nc.gt.nr) then
            call icopy(nr,ir,1,qq(kvt),1)
            call icopy(nc,ic,1,ir,1)
            call icopy(nr,qq(kvt),1,ic,1)
            nn = nc
            nc = nr
            nr = nn
            nn     = iatom1
            iatom1 = iatom2
            iatom2 = nn
            nn  = ib1
            ib1 = ib2
            ib2 = nn
            nn = iorig
            iorig = jorig
            jorig = nn
         end if
         do 60 j=1,nc
            jorb = ic(j)
            call cntrc(qq(ksao),v2(1,jorb),qq(kvt),nbasis,crilow)
            do 50 k=1,nr
               korb = ir(k)
               s12 = ddot(nbasis,qq(kvt),1,v2(1,korb),1)
               qq(ksvd+(j-1)*nr+k-1) = s12
50          continue
60       continue
c.....
c.....   perform a singular value decomposition
c.....
         kdd = igmem_alloc_inf(nc,'vbin.m','hybmak',
     &                         'kdd',IGMEM_DEBUG)
         kud = igmem_alloc_inf(nc*nr,'vbin.m','hybmak',
     &                         'kud',IGMEM_DEBUG)
         kvd = igmem_alloc_inf(nc*nc,'vbin.m','hybmak',
     &                         'kvd',IGMEM_DEBUG)
         kscrd = igmem_alloc_inf(nc,'vbin.m','hybmak',
     &                           'kscrd',IGMEM_DEBUG)
         call svd(nr,nc,qq(ksvd),qq(kdd),.true.,qq(kud)
     &            ,.true.,qq(kvd),ie,qq(kscrd))
         imax = 0
         smax = 0.0d0
         do 70 k=1,nc
c....       take the product of the weights of the original bonding
c....       orbitals and the overlap as criterion
            ip1 = kud + (k-1)*nr + iorig - 1
            ip2 = kvd + (k-1)*nc + jorig - 1
            scrit = dabs(qq(kdd+k-1)) * dabs(qq(ip1)) * dabs(qq(ip2))
            if (scrit.gt.smax) then
c....          acceptable guess
               smax = scrit
               imax = k
            end if
70       continue
         if (imax.eq.0) then
c....       forget it
            write(iwr,11) atoms(iatom1),atoms(iatom2)
11          format(' no hybrid for ',a8,'--',a8,' due to lack of',
     &                                          ' overlap')
         else if (ie.ne.0) then
c....       forget it
            write(iwr,22) atoms(iatom1),atoms(iatom2)
22          format(' no hybrid for ',a8,'--',a8,' due to a sick',
     &                                          ' overlap')
         else
            ibig1 = (imax-1)*nr
            ibig2 = (imax-1)*nc
1            if (ibefo1.ne.0) then
c....          the orbital "resonates" => compromise, that is
c....          add the guesses
               call fmove(v(1,ib1),qq(kvt),nbasis)
            else
               call vclr(qq(kvt),1,nbasis)
            end if
c....       copy the hybrid in the vector array
            do k=1,nr
                call daxpy(nbasis,qq(kud+k+ibig1-1),
     1                     v2(1,ir(k)),1,qq(kvt),1)
            end do
            call fmove(qq(kvt),v(1,ib1),nbasis)
c....       normalise
            nrone=1
            call normt( v(1,ib1),nbasis,nrone,qq(ksao),qq(kvt))
            if (ibefo2.ne.0) then
c....          the orbital "resonates" => compromise, that is
c....          add the guesses
               call fmove(v(1,ib2),qq(kvt),nbasis)
            else
               call vclr(qq(kvt),1,nbasis)
            end if
            do k=1,nc
                call daxpy(nbasis,qq(kvd+k+ibig2-1),
     1                     v2(1,ic(k)),1,qq(kvt),1)
            end do
            call fmove(qq(kvt),v(1,ib2),nbasis)
            nrone=1
            call normt( v(1,ib2),nbasis,nrone,qq(ksao),qq(kvt))
         end if
c
         call gmem_free_inf(kscrd,'vbin.m','hybmak','kscrd')
         call gmem_free_inf(kvd,'vbin.m','hybmak','kvd')
         call gmem_free_inf(kud,'vbin.m','hybmak','kud')
         call gmem_free_inf(kdd,'vbin.m','hybmak','kdd')
         call gmem_free_inf(ksvd,'vbin.m','hybmak','ksvd')
100   continue
      if (iprins.gt.10.or.nscf.eq.0.or.1.eq.1) then
c....    collect the guessed hybrids
         kit = igmem_alloc_inf(2*nbonds*nbasis,'vbin.m','hybmak',
     &                         'kit',IGMEM_DEBUG)
         kittt = igmem_alloc_inf(nbonds,'vbin.m','hybmak',
     &                           'kittt',IGMEM_DEBUG)
         it = 0
         ittt = 0
         do 110 i=1,nbonds
            ib1 = ibond(i,1)
            ib2 = ibond(i,2)
c           if (ib1.lt.0.or.ib2.lt.0) goto 110
            call fmove(v(1,ib1),qq(kit+it*nbasis ),nbasis)
            it = it + 1
            call fmove(v(1,ib2),qq(kit+it*nbasis),nbasis)
            it = it + 1
            ittt = ittt + 1
            call cntrc(qq(ksao),v(1,ib1),qq(kvt),nbasis,crilow)
            qq(kittt-1+ittt) = ddot(nbasis,qq(kvt),1,v(1,ib2),1)
110      continue
         write(iwr,*)(' guessed bonded hybrid pairs :')
         if (iprins.gt.20) call prvc(qq(kit),it,nbasis,dummy,'o','l')
         write(iwr,44) (qq(kittt-1+ij),ij=1,ittt)
44       format(' overlaps of hybrids, per pair :',(t35,13f7.3),/)
         call gmem_free_inf(kittt,'vbin.m','hybmak','kittt')
         call gmem_free_inf(kit,'vbin.m','hybmak','kit')
      end if
c.....make the bond-numbers positive
c     do 120 i=1,nbonds
c        ibond(i,1) = iabs(ibond(i,1))
c        ibond(i,2) = iabs(ibond(i,2))
c120   continue
      call izero(nsingly,isingly,1)
      nsingly = 0
      do i=1,nbonds
         ibond(i,1) = ibond_save(i,1)
         ibond(i,2) = ibond_save(i,2) 
      end do
c
      call gmem_free_inf(kvt,'vbin.m','hybmak','kvt')
      call gmem_free_inf(ksao,'vbin.m','hybmak','ksao')
      return
      end
      subroutine svd(m,n,a,w,matu,u,matv,v,ierr,rv1)
c
c... routine available on SGI. maybe not on others check!
c
      implicit REAL (a-h,o-z), integer(i-n)
c.....
c.....singular value decomposition.
c.....stolen from GNOME
c.....
c
      dimension a(m,n),w(n),u(m,n),v(n,n),rv1(n)
      logical matu,matv
      ierr = 0
      do 100 i=1,m
         do 100 j=1,n
            u(i,j) = a(i,j)
100   continue
      g = 0.0d0
      scale = 0.0d0
      anorm = 0.0d0
      do 300 i=1,n
         l = i + 1
         rv1(i) = scale * g
         g = 0.0d0
         s = 0.0d0
         scale = 0.0d0
         if (i.gt.m) goto 210
         do 120 k=i,m
120      scale = scale + dabs(u(k,i))
         if (scale.eq.0.0d0) goto 210
         do 130 k=i,m
            u(k,i) = u(k,i) / scale
            s = s + u(k,i) * u(k,i)
130      continue
         f = u(i,i)
         g = -dsign(dsqrt(s),f)
         h = f * g - s
         u(i,i) = f - g
         if (i.eq.n) goto 190
         do 150 j=l,n
            s = 0.0d0
            do 140 k=i,m
140         s = s + u(k,i) * u(k,j)
            f = s / h
            do 150 k=i,m
               u(k,j) = u(k,j) + f * u(k,i)
150      continue
190      do 200 k=i,m
200      u(k,i) = scale * u(k,i)
210      w(i) = scale * g
         g = 0.0d0
         s = 0.0d0
         scale = 0.0d0
         if (i.gt.m.or.i.eq.n) goto 290
         do 220 k=l,n
220      scale = scale + dabs(u(i,k))
         if (scale.eq.0.0d0) goto 290
         do 230 k=l,n
            u(i,k) = u(i,k) / scale
            s = s + u(i,k) * u(i,k)
230      continue
         f = u(i,l)
         g = -dsign(dsqrt(s),f)
         h = f * g - s
         u(i,l) = f - g
         do 240 k=l,n
240      rv1(k) = u(i,k) / h
         if (i.eq.m) goto 270
         do 260 j=l,m
            s = 0.0d0
            do 250 k=l,n
250         s = s + u(j,k) * u(i,k)
            do 260 k=l,n
               u(j,k) = u(j,k) + s * rv1(k)
260      continue
270      do 280 k=l,n
280      u(i,k) = scale * u(i,k)
290      anorm = dmax1(anorm,dabs(w(i))+dabs(rv1(i)))
300   continue
      if (.not.matv) goto 410
      do 400 ii=1,n
         i = n + 1 - ii
         if (i.eq.n) goto 390
         if (g.eq.0.0d0) goto 360
         do 320 j=l,n
320      v(j,i) = ( u(i,j) / u(i,l) ) / g
         do 350 j=l,n
            s = 0.0d0
            do 340 k=l,n
340         s = s + u(i,k) * v(k,j)
            do 350 k=l,n
               v(k,j) = v(k,j) + s * v(k,i)
350      continue
360      do 380 j=l,n
            v(i,j) = 0.0d0
            v(j,i) = 0.0d0
380      continue
390      v(i,i) = 1.0d0
         g = rv1(i)
         l = i
400   continue
410   if (.not.matu) goto 510
      mn = n
      if (m.lt.n) mn = m
      do 500 ii=1,mn
         i = mn + 1 - ii
         l = i + 1
         g = w(i)
         if (i.eq.n) goto 430
         do 420 j=l,n
420      u(i,j) = 0.0d0
430      if (g.eq.0.0d0) goto 475
         if (i.eq.mn) goto 460
         do 450 j=l,n
            s = 0.0d0
            do 440 k=l,m
440         s = s + u(k,i) * u(k,j)
            f = (s / u(i,i)) / g
            do 450 k=i,m
               u(k,j) = u(k,j) + f * u(k,i)
450      continue
460      do 470 j=i,m
470      u(j,i) = u(j,i) / g
         goto 490
475      do 480 j=i,m
480      u(j,i) = 0.0d0
490      u(i,i) = u(i,i) + 1.0d0
500   continue
510   do 700 kk=1,n
         k1 = n - kk
         k = k1 + 1
         its = 0
520      do 530 ll=1,k
            l1 = k - ll
            l = l1 + 1
c....       this looks silly, doesn't it  ?
            if (dabs(rv1(l))+anorm.eq.anorm) goto 565
            if (dabs(w(l1))+anorm.eq.anorm) goto 540
530      continue
540      c = 0.0d0
         s = 1.0d0
         do 560 i=l,k
            f = s * rv1(i)
            rv1(i) = c * rv1(i)
            if (dabs(f)+anorm.eq.anorm) goto 565
            g = w(i)
            h = dsqrt(f*f+g*g)
            w(i) = h
            c = g/h
            s = -f/h
            if (.not.matu) goto 560
            do 550 j=1,m
               y = u(j,l1)
               z = u(j,i)
               u(j,l1) = y*c + z*s
               u(j,i) = -y*s + z*c
550         continue
560      continue
565      z = w(k)
         if (l.eq.k) goto 650
         if (its.eq.30) goto 1000
         its = its + 1
         x = w(l)
         y = w(k1)
         g = rv1(k1)
         h = rv1(k)
         f = ((y-z) * (y+z) + (g-h) * (g+h)) / (2.0d0 * h * y)
         g = dsqrt(f*f+1.0d0)
         f = ((x-z) * (x+z) + h*(y/(f+dsign(g,f))-h)) / x
         c = 1.0d0
         s = 1.0d0
         do 600 i1=l,k1
            i = i1 + 1
            g = rv1(i)
            y = w(i)
            h = s * g
            g = c * g
            z = dsqrt(f*f+h*h)
            rv1(i1) = z
            c = f/z
            s = h/z
            f = x*c + g*s
            g = -x*s + g*c
            h = y*s
            y = y*c
            if (.not.matv) goto 575
            do 570 j=1,n
               x = v(j,i1)
               z = v(j,i)
               v(j,i1) = x*c + z*s
               v(j,i) = -x*s + z*c
570         continue
575         z = dsqrt(f*f+h*h)
            w(i1) = z
            if (z.eq.0.0d0) goto 580
            c = f/z
            s = h/z
580         f = c*g + s*y
            x = -s*g + c*y
            if (.not.matu) goto 600
            do 590 j=1,m
               y = u(j,i1)
               z = u(j,i)
               u(j,i1) = y*c + z*s
               u(j,i) = -y*s + z*c
590         continue
600      continue
         rv1(l) = 0.0d0
         rv1(k) = f
         w(k) = x
         goto 520
650      if (z.ge.0.0d0) goto 700
         w(k) = -z
         if (.not.matv) goto 700
         do 690 j=1,n
690      v(j,k) = -v(j,k)
700   continue
      goto 1001
1000  ierr = k
1001  return
      end
      subroutine inform  (nelec,nalpha,ndets,nstruc,pacdet,ndetps,idetps
     &                ,coeff,ncoeff,igroup,nconf,iorth,north,iq,lword)
      implicit REAL  (a-h,o-z) ,  integer   (i-n)
INCLUDE(common/scftvb)
INCLUDE(../m4/common/sizes)
INCLUDE(../m4/common/iofile)
INCLUDE(common/turtleparam)
INCLUDE(common/c8_16vb)
INCLUDE(common/brill)
INCLUDE(common/twice)
INCLUDE(common/hsmattype)
c***********************************************************************
c
c     'inform' reads the information that is written on tape25, by the
c     datatape-program 'crestr'.
c      reorder determinants using symmetry unfo in iorth(input)
c     iq is scratch, lword long
c
      dimension pacdet ( * ),
     &          ndetps ( * ),
     &          idetps ( * ),
     &          coeff  ( * ),
     &          igroup ( 5 , * ),
     &          iorth  ( * ),
     &          iq(lword)
c
c                          read variables
c
_IF(atmol)
      rewind 25
      read(25) nelec,nalpha,nconf,ndets,nstruc,maxdet,maxstr,
     &         nwpack,ncoeff,imax
_ELSE
      call wr15vb(nelec,nalpha,nconf,ndets,nstruc,maxdet,maxstr,
     &            nwpack,ncoeff,imax,ndum,ndoub,'read')
_ENDIF
      kscr = 1
      kidet = kscr + maxdet
      kdoub = kidet + maxdet * nelec
      klast = kdoub + ndoub
      kscrgr = klast
      klast = klast + nconf
      if(klast.ge.lword)call corfait(klast,lword,'inform')
c      print *,' nelec,nalpha,nconf,ndets,nstruc,maxdet,maxstr,'
c     &       ,' nwpack,ncoeff,imax'
c      print *, nelec,nalpha,nconf,ndets,nstruc,maxdet,maxstr,
c     &         nwpack,ncoeff,imax
c
c               read number of determinants per structure
c
_IF(atmol)
      read (25) (ndetps(i),i=1,nstruc)
_ELSE
      call readis(ndetps,nstruc,num8)
_ENDIF
c      print*,'ndetps :',(ndetps(i),i=1,nstruc)
c
c                      read determinant-numbers
c
      ndettot = 0
      do 10 i=1,nstruc
         ndettot = ndettot + ndetps(i)
10    continue
_IF(atmol)
      read(25) (idetps(i),i=1,ndettot)
_ELSE
      call readis(idetps,ndettot,num8)
_ENDIF
c      print*,'idetps :', (idetps(i),i=1,ndettot)
c
c                         read coefficients
c
_IF(atmol)
      read(25) (coeff(i),i=1,ndettot)
_ELSE
      call reads(coeff,ndettot,num8)
_ENDIF
c      print*,'coeff :', (coeff(i),i=1,ndettot)
c
c                read number of determinants per configuration
c
_IF(atmol)
      read(25) (igroup(1,i),i=1,nconf)
_ELSE
      call readis(iq(kscrgr),nconf,num8)
      do i=1,nconf
         igroup(1,i) = iq(kscrgr+i-1)
      end do
_ENDIF
c      print*,'igroup 1 :', (igroup(1,i),i=1,nconf)
c
c                 read number of structures per configuration
c
_IF(atmol)
      read(25) (igroup(2,i),i=1,nconf)
_ELSE
      call readis(iq(kscrgr),nconf,num8)
      do i=1,nconf
         igroup(2,i) = iq(kscrgr+i-1)
      end do
_ENDIF
c      print*,'igroup 2 :', (igroup(2,i),i=1,nconf)
c
      it = 1
      do 62 i=1,nconf
         igroup(3,i) = 0
         do 61 j=1,igroup(2,i)
            igroup(3,i) = igroup(3,i) + ndetps(it)
            it = it + 1
61       continue
62    continue
c      print*,'igroup 3 :', (igroup(3,i),i=1,nconf)
c
c              read packed slater determinants
c
      nwords = 0
      do 66 i=1,nconf
         nwords = nwords + ((igroup(1,i) * nelec - 1)/(64/n8_16)) + 1
66    continue
c Initialisation of brillioun state optimisation type (default Super CI)
c 1 - Super CI
c 2 - Perturbation theory (aNR) - thesis J. Verbeek
c 3 - Fock matrix (1st column and diagonal)
      do 67 i=1,nconf
         igroup(4,i) = 1
67    continue
_IF(atmol)
      read(25) (pacdet(i),i=1,nwords)
_ELSE
      call reads(pacdet,nwords,num8)
_ENDIF

      if (hsmatt.eq.'vbfunction') then
c.....
         k7doub =  kscra7vb('k7doub',ndoub,'i','r')
         call readi(iq(kdoub),ndoub,k7doub,num8)
         it  = 1
         it2 = 1
         nnelec=nelec-ndoub*2
         do i=1,nconf
            call izero(nelec*igroup(1,i),iq(kidet),1)
            call unpack(pacdet(it),n8_16,iq(kidet),nelec*igroup(1,i))

c      do ij=1,igroup(1,i)
c        print *,ij,' det: ',(iq(kidet+iij-1+ij*nelec-nelec),
c     &          iij=1,nelec*igroup(1,i))
c      end do
            call deldoub(iq(kdoub),iq(kidet),nelec,nalpha,
     &                  ndoub,igroup(1,i))
            call pack(pacdet(it2),n8_16,iq(kidet),nnelec*igroup(1,i))
            it  = it + ((igroup(1,i)*nelec-1)/(64/n8_16))+1
            it2 = it2 + ((igroup(1,i)*nnelec-1)/(64/n8_16))+1
         enddo
         nelec = nnelec
         nalpha = nalpha - ndoub
         nwpack = ((igroup(1,1)*nnelec-1)/(64/n8_16))+1
         imax = imax - ndoub
      endif
c.....
      it  = 1
      it2 = 1
c.....
c.....reorder the determinants so that orbitals occur in orthogonality
c.....groups
c.....
      if (nitscf.eq.1.and.iprins.gt.1000) then
         write(iwr,999)
999      format(/,' the determinants and their spin-coefficients :')
      end if
      do 444 i=1,nconf
         call izero(nelec*igroup(1,i),iq(kidet),1)
         call unpack(pacdet(it),n8_16,iq(kidet),nelec*igroup(1,i))
         ntot = igroup(3,i)
c         print *,'i: ',i,' group 1: ',igroup(1,i)
c         print *,'group 3: ',igroup(3,i)
         call reorb(iq(kidet),igroup(1,i),nelec,nalpha,iorth,iq(kscr),
     &                               coeff(it2), idetps(it2),ntot)
         call pack(pacdet(it),n8_16,iq(kidet),nelec*igroup(1,i))
         if (nitscf.eq.1.and.iprins.gt.50000) then
            write(iwr,234)i
234         format(/,' configuration number',i3)
            ipr = kidet
            do 987 j=1,igroup(1,i)
             write(iwr,645) coeff(it2+j-1),(iq(ij),ij=ipr,ipr+nalpha-1),
     &                             (iq(ij),ij=ipr+nalpha,ipr+nelec-1)
645          format(1x,f7.4,(40i3),' - ',(40i3))
             ipr = ipr + nelec
987         continue
         end if
         it  = it + ((igroup(1,i)*nelec-1)/(64/n8_16))+1
         it2 = it2 + ntot
444   continue
      return
      end
      subroutine lowdins(c,sao,s,cmsq,lambda,vecmt,n,m,epsss)
c
c...  begin of routines, stolen from servec for lowdin orthogonalisation
c
      implicit REAL  (a-h,o-z) , integer   (i-n)
c
c
c     subroutine lowdin performs a lowdin (symmetric) orthonormalization
c     of a given set of n vectors of dimension m
c     input are the vector-array c and the triangular primitive
c     overlap-matrix sao
c     the results are returned in the vector-input array
c
c     arrays
c            c      : n*m
c            sao    : m*(m+1)/2
c            s      : n*(n+1)/2
c            cmsq   : n*n
c            lambda : n
c            vecmt  : m
c
c
      dimension c(1),sao(1),s(1),lambda(1),cmsq(1),vecmt(1)
      integer slen
      REAL lambda
c
c     form the overlap-matrix s
c
      call fmos(s, sao, c, vecmt, n, m, epsss)
      k = n*(n+1)/2
c
c     diagonalize the s-matrix
c
      call jaclow(s,cmsq,n,k,epsss)
c
c     form s**(-.5)
c
c     make a row of the negative squareroot of the eigenvalues
c
      k = 0
      do 70 i=1,n
         k = k + i
         if (dabs(s(k)).lt.epsss) then
c....       eigenvector effectively zero
            lambda(i) = 0.0d0
         else
            lambda(i) = 1/dsqrt(s(k))
         end if
70    continue
c
c     perform the matrix-multiplication
c
      k = 0
      do 90 i=1,n
      kcmin = i - n
      do 80 j=1,n
      kcmin = kcmin + n
   80 vecmt(j) = cmsq(kcmin)*lambda(j)
      do 90 j=1,i
      kl = j - n
      res = 0.0d0
      do 85 ij =1,n
      kl = kl + n
   85 res = res + vecmt(ij)*cmsq(kl)
      k = k + 1
   90 s(k) = res
c
c     bring the lowdin-transformation-matrix on rectangle-form
c
      k = 0
      l = 0
      iend = n - 1
      do 110 i=1,iend
      do 100 j=1,i
      k = k + 1
      l = l + 1
  100 cmsq(l) = s(k)
      jbeg = i + 1
      kl = k + i
      do 110 j=jbeg,n
      l = l + 1
      cmsq(l) = s(kl)
  110 kl = kl + j
      do 120 i=1,n
      k = k + 1
      l = l + 1
  120 cmsq(l) = s(k)
c
c     transformation of the vectors
c
      call matmpl(c,cmsq,vecmt,n,m)
c
      return
      end
c***********************************************************************
      subroutine dumporb(v,m,n,ndim)

      implicit REAL  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
c
INCLUDE(../m4/common/sizes)
INCLUDE(../m4/common/prints)
INCLUDE(../m4/common/iofile)
INCLUDE(../m4/common/runlab)
INCLUDE(common/dumpvbmo)
c
      dimension v(ndim,*)
c
c...  Routine for dumping VB MOs to separate file, in order with
c...  MOLDEN printing format
c
      max = 10
      if (oprint(20)) max=7
      imax = 0
  100 imin = imax+1
      imax = imax+max
      if (imax .gt. m) imax = m
        write(vbmounit,1000)
        write(vbmounit,8028) (i-i,i = imin,imax)
        write(vbmounit,1000)
        write(vbmounit,8028) (i,i = imin,imax)
        write(vbmounit,1000)
        do 120 j = 1,n
          write(vbmounit,8048) j,zbflab(j),(v(j,i),i = imin, imax)
 120    continue
      if (imax .lt. m) go to 100
      return

 8028 format(17x,10(3x,i3,3x))
 8048 format(i5,2x,a10,10f9.4)
 1000 format(/)

      end
c***********************************************************************
      subroutine prvc (c,m,n,eps,ch,test)
c
      implicit REAL  (a-h,o-z) , integer   (i-n)
_IFN(atmol)
INCLUDE(../m4/common/iofile)
_ENDIF
c
c   routine for printing of column-vectors and eigenvalues.
c   now heavily changed with respect to original (esp. parameters)
c   # columns fixed to 10 // 6 decimal places
c
c            c  : matrix of coefficients
c            m  : x columns to be printed
c            n  : dimension of c(n,m) and eps(n)
c            eps: vector of eigenvalues/occupations
c            ch : 'value' or 'v' => values; 'order' or 'o' => columns 
c           test: 'label' ior 'l' => gamess labels if possible
c
c
      dimension c(n,m), eps(m) 
      character*(*) test,ch
      dimension sub(2)
      character*5 sub
      data sub /'-----','-----'/
c     

      if (ch.eq.'value'.or.ch.eq.'v') then
         if (test.eq.'label'.or.test.eq.'l') then
_IFN(atmol)
            call prev(c,eps,m,n,n)
         else  
_ENDIF
            ich = 1
            go to 1
         end if
      else if (ch.eq.'order'.or.ch.eq.'o') then 
         if (test.eq.'label'.or.test.eq.'l') then
_IFN(atmol)
            call prsql(c,m,n,n)
         else  
_ENDIF
            ich = 0
            go to 1
         end if
      else 
         call vberr(' wrong prvc call')
      endif 

      return
c
1     ncc = 10
      nbl = (m-1)/ncc
      nlast = m-ncc*nbl
      nbl = nbl+1
      nc1 = 0
      do 60 ib=1,nbl
      if ( ib .eq. nbl ) ncc=nlast
      nc0 = nc1+1
      nc1 = nc1+ncc
      if (ich.ne.0) write(iwr,10) ( eps(ic), ic=nc0,nc1 )
      if (ich.eq.0) write(iwr,11) ( ic, ic=nc0,nc1 )
   10 format(/,7x,10f12.6)
   11 format(/,7x,10(3x,i4,3x))
      write(iwr,20) ( sub, i=1,ncc )
   20 format(8x,10(1x,a5,a5,1x))
      write(iwr,30)
   30 format(1h )
      do 50 ia=1,n
      write(iwr,40) ia, ( c(ia,ic), ic=nc0,nc1 )
   40 format(2x,i3,2x,10f12.6)
   50 continue
   60 continue
      write(iwr,30)
c
      return
      end

c***********************************************************************
      subroutine prvcacc(c,m,n,isect)
      implicit REAL  (a-h,o-z) , integer   (i-n)
      dimension c(m,n)
INCLUDE(common/scftvb)

      ncc = 4
      nbl = (m-1)/ncc
      nlast = m-ncc*nbl
      nbl = nbl+1
      nc1 = 0

      do ib=1,nbl
      if ( ib .eq. nbl ) ncc=nlast
      nc0 = nc1+1
      nc1 = nc1+ncc

      do ia=1,n
        write(isect,'(I3,a,4F22.16)') ia,'    ',(c(ia,ic), ic=nc0,nc1)
      end do

      write(isect,*)
      end do

      return
      end

c***********************************************************************

      subroutine rdist1(istr,nn)
c
      implicit REAL  (a-h,o-z) , integer   (i-n)
c
c...  read a string numbers (including 'to') on 1 card
c
      dimension istr(*)
c
      logical to
c
      character*4 iend,ito,ispace,itest
      data iend,ito,ispace/'end ','to','    '/
c
      nn = 0
      to = .false.
c
10    if (to) then
         call inpi(ne)
         if (ne.le.0) call vberr('  missing 2nd arg in ..to.. string')
         if (nn.le.0) call vberr('  missing 1st arg in ..to.. string')
         nb = istr(nn)
         nn = nn - 1
         do 20 ii=nb,ne
            nn = nn + 1
            istr(nn) = ii
20       continue
         to = .false.
      else
         call inpa4(itest)
         if (itest.eq.iend.or.itest.eq.ispace) then
            return
         else if (itest.eq.ito) then
            to = .true.
         else
            call cinput(jrec,jump,-1)
            jrec = jrec - 1
            nn = nn + 1
            call inpi (istr(nn))
         end if
      end if
c
      goto 10
c
      return
      end
      subroutine rdistr(istr,nn,max)
c
      implicit REAL (a-h,o-z), integer (i-n)
c
c...  read a string numbers (including 'to') until end is encountered
c...  max is maximum # numbers to be read
c
      character*8 itest,iend,ito
      dimension istr(*)
      logical to
      data iend,ito/'end','to'/
c
      nn = 0
      to = .false.
       do i=1,max
          istr(i) = 0
       end do
c
10    call rchar(itest)
      if (itest.eq.iend) go to 100
      if (to) go to 20
      if (itest.ne.ito) go to 20
      to = .true.
      go to 10
c
20    call cinput(jrec,jump,-1)
      call inpi(ne)
      if (to) go to 30
      nn = nn + 1
      istr(nn) = ne
      go to 10
c
30    nb = istr(nn)
      nn = nn - 1
      do 40 ii=nb,ne
         nn = nn + 1
40    istr(nn) = ii
c
      to = .false.
      go to 10
c
100   if (to) nn = max
      if  (nn.gt.max) call vberr('too many numbers in rdistr')
c
      return
      end
      subroutine rdistr_plus(istr,nstr,rstr,maxstr,maxnn,nn)
c
      implicit REAL (a-h,o-z), integer (i-n)
c
c...  read a string numbers (including 'to') until end is encountered
c...  this is a special version for reading combinations (of vectors)
c...  now input like 0.3484 * 4 + 0.3 * 7   8 9 is allowed
c..
c...  'to' 'end' is disabled (to unclear anyway) 
c...  
c...   returns :
c...     istr : integers read
c...     nstr : # of (vectors) to be added including first one
c...     rstr : coefficients of (vectors) 
c...     nn   : number of (vectors) to be treated
c
      character*8 test
      dimension rstr(maxstr,maxnn),istr(maxstr,maxnn),nstr(maxnn)
c
      nn = 0
      iplus = 0
c
10    call rchar(test)
c
      if (test.eq.'end') then
         return
c
      else if (test.eq.'to') then
         if (nstr(nn).gt.1) call caserr(' no to with + ')
         nb = istr(1,nn)
         call inpi(ne)
         nn = nn - 1
         do ii=nb,ne
            nn = nn + 1
            nstr(nn) = 1
            istr(1,nn) = ii
            nstr(nn) = 1
            rstr(1,nn) = 1.0d0
         end do
c
      else if (index(test,'.').ne.0) then
         if (iplus.eq.0) then
            nn = nn + 1
            nstr(nn) = 0
            if (nn.gt.maxnn) 
     1         call caserr('istr overflow in rdistr_plus')
         end if
c...     floating point number (mult factor) * integer
         nstr(nn) = nstr(nn) + 1
         if (nstr(nn).gt.maxstr) 
     1      call caserr('rstr overflow in rdistr_plus')
         call cinput(jrec,jump,-1)
         call inpf(rstr(nstr(nn),nn))
         if (iplus.eq.-1) rstr(nstr(nn),nn) = rstr(nstr(nn),nn)*(-1.0d0)
         call inpa4(test)
         if (test.ne.'*') call caserr('mult-factor not followed by *')
         call inpi(istr(nstr(nn),nn))
         iplus = 0
      else if (test.eq.'+') then
         iplus = +1
      else if (test.eq.'-') then
         iplus = -1
      else
         if (iplus.eq.0) then
            nn = nn + 1
            nstr(nn) = 0
            if (nn.gt.maxnn) 
     1         call caserr('istr overflow in rdistr_plus')
         end if
         call cinput(jrec,jump,-1)
         nstr(nn) = nstr(nn) + 1
         if (nstr(nn).gt.maxstr) 
     1      call caserr('rstr overflow in rdistr_plus')
         call inpi(istr(nstr(nn),nn))
         rstr(nstr(nn),nn) = 1.0d0
         if (iplus.eq.-1) rstr(nstr(nn),nn) = -1.0d0
         iplus = 0
      end if
c
      go to 10
c

      return
      end
c***********************************************************************
      subroutine restror(v,nvec,nbasis,qq)
c
      implicit REAL (a-h,o-z), integer (i-n)
c
      dimension v(nbasis,nvec),qq(*)
c
INCLUDE(../m4/common/sizes)
INCLUDE(../m4/common/gmempara)
INCLUDE(common/turtleparam)
      common /forbex/ nforbi,iforbi(maxact*maxact,2),
     &                nionj,ionj(maxact,2),
     &                nsetor,isetor(maxact+1,maxact),
     &                nionres,ionrest(maxact)
c
INCLUDE(common/infato)
INCLUDE(common/vbcri)
c.....
c.....perform orthogonalisations as defined in forbex
c.....
      ks = igmem_alloc_inf(nbasis*(nbasis+1)/2,'vbin.m','restror','ks',
     &                     IGMEM_DEBUG)
      call get1e(qq(ks),dummy,'s',qq(ks))
c.....
c.....i on j
c.....
      kw = igmem_alloc_inf(nbasis,'vbin.m','restror','kw',IGMEM_DEBUG)
      do 10 ij=1,nionj
         i = ionj(ij,1)
         j = ionj(ij,2)
         call schmids(v(1,i),v(1,j),qq(ks),qq(kw),
     &                1,1,nbasis,cridep,'norm')
10    continue
c.....
c.....i on rest ; exclude i of course
c.....
      do 20 i=1,nionres
         nvecc = ionrest(i) - 1
         call schmids(v(1,ionrest(i)),v,qq(ks),qq(kw),
     &                1,nvecc,nbasis,cridep,'norm')
         nvecc = nvec - ionrest(i) + 1
         call schmids(v(1,ionrest(i)),v(1,ionrest(i)+1),qq(ks),qq(kw),
     &                1,nvecc,nbasis,cridep,'norm')
20    continue
c.....
c.....sets
c.....
      do 30 i=1,nsetor
         nset = isetor(1,i)
         kvc = igmem_alloc_inf(nbasis*nset,'vbin.m','restror','kvc',
     &                         IGMEM_DEBUG)
         call vecgat(v,qq(kvc),nbasis,isetor(2,i),nset)
         ksmo = igmem_alloc_inf(nset*(nset+1)/2,'vbin.m','restror',
     &                          'ksmo',IMEM_DEBUG)
         kcmsq = igmem_alloc_inf(nset*nset,'vbin.m','restror',
     &                           'kcmsq',IMEM_DEBUG)
         klambda = igmem_alloc_inf(nset,'vbin.m','restror',
     &                             'klambda',IMEM_DEBUG)
         kvecmt = igmem_alloc_inf(nbasis,'vbin.m','restror',
     &                            'kvecmt',IMEM_DEBUG)
c
         call lowdins(qq(kvc),qq(ks),qq(ksmo),qq(kcmsq),qq(klambda),
     &                qq(kvecmt),nset,nbasis,crilow)
         call vecsca(v,qq(kvc),nbasis,isetor(2,i),nset)
c
         call gmem_free_inf(kvecmt,'vbin.m','restror','kvecmt')
         call gmem_free_inf(klambda,'vbin.m','restror','klambda')
         call gmem_free_inf(kcmsq,'vbin.m','restror','kcmsq')
         call gmem_free_inf(ksmo,'vbin.m','restror','ksmo')
         call gmem_free_inf(kvc,'vbin.m','restror','kvc')
30    continue
c
      call gmem_free_inf(kw,'vbin.m','restror','kw')
      call gmem_free_inf(ks,'vbin.m','restror','ks')
c
      return
      end

c***********************************************************************
      subroutine vecman(v,nbasis,nvec,not)
      implicit REAL (a-h,o-z), integer(i-n)
c.....
c.....reads vectors from file, not is the number of dummy columns
c.....before the actual vectors (numbering etc.).
c.....
INCLUDE(../m4/common/iofile)
      dimension v(nbasis,nvec)
      character*80 filen
      character*1 line(132)
      character*132 line2
      equivalence (line,line2)
      logical space
      nmo = 0
      if (nvec.lt.1) then
         write(iwr,11)
11       format(/,' please specify the number of vectors that are to',
     &            ' be read from file')
         call vberr('     try again')
      end if
10    continue
      ifile = ird  
c     if (nmo+ncurr.gt.nvec) ncurr = nvec - nmo
20    do 30 i=1,nbasis
123      read(ifile,33,end = 1000) line
33       format(132a1)
         if (line(1).eq.'*'.or.line(1).eq.'?') goto 123
         space = .true.
         ndumm = 0
         do 25 j=1,132
            if (line(j).ne.' ') then
               if (space) then
                  ndumm = ndumm + 1
                  space = .false.
               end if
               if (ndumm.le.not) line(j)=' '
            else
               space = .true.
            end if
25       continue
c....    if we get here it was a blank line
         ncurr = ndumm - not
         if (ndumm.eq.0) goto 123
         read(line2,*,end=1002,err=1002) (v(i,nmo+j),j=1,ncurr)
30    continue
      nmo = nmo + ncurr
      if (nmo.eq.nvec) then
        return
      end if
      goto 20
1000  write(iwr,1001) i,nmo+1,ncurr,line
1001  format(' ** end of file or error** on vectors manual ',/
     1       ' at column ',i5,' and first mo',i5,' # mos ',i5,/,
     2       ' last line  ',/,1x,132a1)
      call vberr(' eof on vecman')
1002  write(iwr,1001) i,nmo+1,ncurr,line
      call vberr(' eof on internal read in vecman ')
      return
      end
      subroutine veccomb(q,qq,qqt,map,map2,nstr,rstr,istr,
     1                   ndim,ncol,maxcol,maxstr)
c
c...   set up start-vectors from combined sectors or input sets
c...   options are added; routine is quickly getting out of hand
c
      implicit REAL  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
c
      dimension q(ndim,maxcol),qq(ndim,maxcol),qqt(ndim,maxcol) 
      dimension map(maxcol),map2(maxcol)
      dimension nstr(maxcol),rstr(maxstr,maxcol),istr(maxstr,maxcol)
INCLUDE(common/scftvb)
INCLUDE(common/vbcri)
INCLUDE(common/turtleparam)
INCLUDE(common/aivb)
INCLUDE(../m4/common/iofile)
INCLUDE(../m4/common/sector)
      tresh = crihyb*100.0d0
c
      numdu_save = 0
      ibldu_save = 0
c
      nv = 0
1     call input
      call inpa4(ytest)
      if (ytest.eq.'sect'.or.ytest.eq.'file') then
         if (ytest.eq.'file') then
            call inpa4(ytest)
            call inpi(iblk)
            call filec(ytest,iblk,iunit,irep)
            if (irep.ne.0) call caserr('wrong dumpfile in veccomb')
            numdu_save = numdu
            ibldu_save = iblkdu
            call revind
            call secini(iblk,iunit)
         end if 
         call inpi(isec)
         call rdistr_plus(istr,nstr,rstr,maxstr,maxcol,nn)
         call getqvb(qq,nbasis,ncol,isec,'print')
c
         if (nbasis.ne.ndim) then
            call input
            call inpa4(ytest)
            if (ytest.ne.'extr')call caserr('no expected extra in comb')
            call rdistr(map2,nne,maxcol)
            call vbextr(qq,ndim,qqt,nbasis,map2,nne)
            nbasis = nbasis + nne
            if (nbasis.ne.ndim) call caserr('veccom nbasis ne ndim')
         end if
         if (nv+nn.gt.maxcol) call caserr('# vectors large (veccomb)')
         do i=1,nn
            nv = nv + 1
            call vclr(q(1,nv),1,nbasis)
            do j=1,nstr(i)
                 call daxpy(nbasis,rstr(j,i),qq(1,istr(j,i)),1,
     1                      q(1,nv),1)
            end do
         end do
         if (numdu_save.ne.0) then
            call secini(ibldu_save,numdu_save)
            numdu_save = 0
            ibldu_save = 0
         end if
         write(iwr,603) nv
603      format(1x,'*** now ',i6,' vectors ***')
      else if (ytest.eq.'manu') then
         call inpi(nvec)
         call inpi(not)
         if (nv+nvec.gt.maxcol) call caserr('# vectors large (veccomb)')
         call vecman(q(1,nv+1),nbasis,nvec,not)
         nv = nv + nvec
         write(iwr,603) nv
      else if (ytest.eq.'dump') then
c.....   dump the vectors now read on the dumpfile 
         call inpi(isecv)
         write(iwr,605) nv,isecv
605      format(' *** dump current vectors (#',i4,') to sect',i3,' ***')
         call putqvb(q,nbasis,nv)
         isecv = 0
      else if (ytest.eq.'end') then
         call inpa4(ytest)
         if (ytest.ne.'comb'.and.ytest(1:1).ne.' ')
     1    call caserr('unrecognised end in vectors comb')
         ncol = nv
         return
      else if (ytest.eq.'perm'.or.ytest.eq.'mope') then
c...     permute mo's now in core
         call izero(nv,map,1)
         call rdistr(map,nn,maxcol)
         call coperm(map,qq,nv,iwr)  
         write(iwr,606) (map(i),i=1,nv)
606      format(' *** Permutation of current vectors ',(t37,20i4))
         do i=1,nv
            if (map(i).ne.i) then
c...   save vector on position i
               call dcopy(ndim,q(1,i),1,qq,1)
               k = i
c...   copy vector on position map(k) to position k
20             call dcopy(ndim,q(1,map(k)),1,q(1,k),1)
               l = map(k)
               map(k) = k
               if (map(l).ne.i) then
                  k = l 
                  go to 20
               end if
c...  copy original vector on position i to podition just freed
               call dcopy(ndim,qq,1,q(1,l),1)
               map(l) = l
            end if
         end do
      else if (ytest.eq.'elim') then
c...     eliminate vectors currently in core having cofficients at specified places
         call rdistr(map,nn,maxcol)
         write(iwr,601) map(1),(map(i),i=2,nn)
601      format(' eliminate exactly',i6,' vectors nonzero at aos',
     1          (t47,10i4))
         nnv = 0
         ib = 1
         qqq = 100.0d0
30       do i=ib,nv
           do k=2,nn
             if (dabs(q(map(k),i)).gt.tresh) then
               qqq = dmin1(dabs(q(map(k),i)),qqq)
c...           eliminate i
               nnv = nnv + 1
               do j=i+1,nv
                 call dcopy(ndim,q(1,j),1,q(1,j-1),1)
               end do
               ib = i
               nv = nv - 1
               go to 30
             end if
           end do
         end do
         write(iwr,602) nv,qqq
602      format(1x,'***',i6,' vectors remain',
     1          ' ,smallest eliminate ',1pd10.2,'***')
         if (nnv.ne.map(1))call caserr(' wrong estimate of elimination')
      else if (ytest.eq.'scre'.or.ytest.eq.'clea') then
c...     screen vectors (i.e. eliminate small elements)
         call inpf(qxx)
         call inpi(i)
         qxx = qxx * 10**(i*1.0d0)
         qqq = 0.0d0
         do i=1,nv
          do j=1,ndim
             if (dabs(q(j,i)).lt.qxx) then
                qqq = dmax1(qqq,dabs(q(j,i)))
                q(j,i) = 0.0d0
             end if
          end do
         end do
         write(iwr,604) qxx,qqq
604      format(' *** clear all elements under ',1pd12.5,' in set',
     1          ', max was ',1pd12.5,'  ***')
      else
         call caserr('illegal entry in vectors comb')
      end if
c
      go to 1
c
      end
      subroutine vbextr(qf,ndim,qi,nbasis,map,nn)
c
c...  simple  extra routine (qi = auxilary)
c...  nbasis < ndim
c     only produces nbasis vectors 
c
      implicit REAL (a-h,o-z), integer (i-n)
      dimension qf(ndim,*),qi(nbasis,*),map(nn)
c
      call dcopy(nbasis*nbasis,qf,1,qi,1)
      call vclr(qf,1,ndim*ndim)
c
      do i=1,nbasis
         k = 1
         do j=1,ndim
            if (map(k).eq.j) then
               k = k + 1
            else
               qf(j,i) = qi(j-k+1,i)
            end if
         end do
      end do
c
      return
      end
      subroutine showstat
c
      implicit REAL (a-h,o-z), integer (i-n)
      parameter (n35=35)
      common /stats_vb/ noccurr(n35)
INCLUDE(../m4/common/iofile)
c      
_IF(parallel)
      call pg_igop(3002,noccurr,n35,'+')
_ENDIF
      write(iwr,'(/,A20)') 'statistics of cases '
      write(iwr,'(A20)') '--------------------'
      write(iwr,'(A20)') 'case     occ.  perc.'
      nocctot=0
      nocctwo=0
      do 43 i=1,n35
        nocctot = nocctot + noccurr(i)
        if (i.gt.6) nocctwo=nocctwo+noccurr(i)
43    continue
      noccone=nocctot-nocctwo-noccurr(1)
      do 44 i=1,n35 
        write(iwr,'(I3,I10,F7.2)') i,noccurr(i),1.0D2*noccurr(i)/nocctot
44    continue  
      write(iwr,'(A20)') '--------------------'
      write(iwr,'(A8,I10,F7.2)') 'nsing 0:',
     &       noccurr(1),1.0D2*noccurr(1)/nocctot
      write(iwr,'(A8,I10,F7.2)') 'nsing 1:',
     &       noccone,1.0D2*noccone/nocctot
      write(iwr,'(A8,I10,F7.2)') 'nsing 2:',
     &       nocctwo,1.0D2*nocctwo/nocctot
c
      return
      end
      subroutine cgamtovb(numed3,iblk3)
c
      implicit REAL  (a-h,p-w),integer   (i-n),logical    (o)
c
INCLUDE(common/tractlt)
c     common/sector/numdu,iblkdu
INCLUDE(../m4/common/sizes)
INCLUDE(../m4/common/sector)
INCLUDE(../m4/common/infoa)
c
      nbasis = num
      iblk3 = iblkdu
      numed3 = numdu
      return
      end

c***********************************************************************
      subroutine vbfirst(v,maxvec)
c
      implicit REAL  (a-h,o-z) ,  integer   (i-n)
c
      common/bypass/ index4,ihmat,idavid
      common/title/ titel,chtime(4)
      character*80 titel
      character*8 chtime,harmvb
      character*1 ch1,ch2
c
c...  hsmat and bsmat
c
INCLUDE(common/hsinfo)
INCLUDE(common/vbcri)
c
c...   davidson
c
      common/davcom/eshift,thresj,thresd,maxcyc,maxdav,maxsel,iprind,
     *              alter,firsts
      common/davctx/mode
      character*4 mode
      logical alter,firsts
c
c...   scf  (initialised in scfin)
c
INCLUDE(common/scftvb)
INCLUDE(../m4/common/sizes)
INCLUDE(../m4/common/iofile)
INCLUDE(../m4/common/vcore)
INCLUDE(common/turtleparam)
INCLUDE(common/twice)
INCLUDE(common/first_vb)
INCLUDE(common/brill)
      common/davscf/eshsst,thrssj,thrssd,maxssc,maxssv,maxssl,iprssd,
     *              alssr,firsss
      common/davssx/modess
      character*4 modess
      logical alssr,firsss
c
c...  for dynamic memory allocation debugging purpose, contains IGMEM_ vars.
INCLUDE(../m4/common/gmemdata)
INCLUDE(../m4/common/gmempara)
c
c...   4-index common-blocks
c
INCLUDE(common/splice)
INCLUDE(common/tractlt)
      common/junke/maxt,ires,ipass,nteff,npass1,npass2,
     *             nbuck,mloww,mhi,ntri,iacc
      common/discc/ ied(16)
INCLUDE(../m4/common/mapper)
INCLUDE(../m4/common/harmon)
INCLUDE(../m4/common/tran)

      character*4 ied
c...  scratch  for mix-directive
      parameter (mxmix=20, mxorg=50)
      common /scrp/ cmix(mxmix,mxorg),imix(mxmix,mxorg),nopomx(mxorg)
c...
c
c...  cases for matrix type analysis
c
      common /cases/ icase(2*ncases),times(ncases),token(3*ncases),
     &               lcases
      logical lcases
c
      logical new2el
      common/new2el/ new2el
c
c...
c...  atonam for names of atoms (for mulliken population analysis)
c...  nopa contains the # orbitals per atom, iopa contains them
c...  mullik indicates if the analysis is requested. iacat contains
c...  scf orbitals per atom, nacat the number of the, hybry indicates
c...  wether the automatic scf-procure should confine itself to
c...  atomic optomisation
c...
INCLUDE(common/infato)
      common/energi/enatom,nstate
c       ... scf etc
      common/gjs   /nirr,mult(8,8),isymao(maxorb),isymmo(maxorb)
     1              ,irr,iss,nstart(8),nfin(8),nmc
c
      logical dump,lschmi,llowdi
      common/comvbin/dump,lschmi,llowdi,nscset,nloset,ntemp,nmix
c
c.....mrsdci
c     jstat control integer for doubles selection:
c
c     jstat < 0 : mrsdci without selection
c     jstat = 0 : no mrsdci (initialised in inivb)
c     jstat = 1 : mrsdci with doubles selection, phase 1
c     jstat = 2 : mrsdci with doubles selection, phase 2
c
      common/logvb/logcor,logact
      logical logcor,logact
c
c...
      dimension v(*)
c...
      parameter (nkey=37)
      character*4 idd(nkey),space,iby(3),iddp,test
c
      kkk = max(nbasis,nsa+ncore+ntemp)
      kkk = max(ncol,kkk)*nbasis
         if (.not.ofirst) go to 1000
      if (lschmi.or.llowdi) then
         maxx = 0
         do 1122 i=1,nscset+nloset
1122     maxx = max(maxx,iex(1,i))

         call get1e(Q(ks),dummy,'s',Q(ks))
         if (lschmi) then
_IF(debug)
            write(iwr,*)  '(1) VBFIRST - LSCHMI'
_ENDIF
            do 4 k=1,nscset
               nv = iex(1,k)
               do 2 i=1,nv
                  call fmove(Q( (iex(i+1,k)-1)*nbasis+1 ),
     &                       Q(kkk+(i-1)*nbasis+1),nbasis)
2              continue
               call normvc(Q(kkk+1),Q(ks),Q(kscr),
     +                     nbasis,nv,cridep)
               do 3 i=1,nv
                  call fmove(Q(kkk+(i-1)*nbasis+1),
     &                       Q( (iex(i+1,k)-1)*nbasis+1 ),nbasis)
3              continue
4           continue
         end if 
         if (llowdi) then
_IF(debug)
            write(iwr,*) '(2) VBFIRST - LLOWDI'
_ENDIF
            do 7 k=1,nloset
               nv = iex(1,nscset+k)
               do 5 i=1,nv
                  call fmove(Q( (iex(i+1,k+nscset)-1)*nbasis+1 ),
     &                       Q(kkk+(i-1)*nbasis+1),nbasis)
5              continue
               ksmo = kscr
               kcmsq = ksmo + nv*(nv+1)/2
               klambda = kcmsq + nv*nv
               kvecmt = klambda + nv
               call lowdins(Q(kkk +1),Q(ks),
     +                      Q(ksmo),Q(kcmsq),
     +                      Q(klambda),Q(kvecmt),
     +                      nv,nbasis,crilow)
               do 6 i=1,nv
                  call fmove(Q(kkk+(i-1)*nbasis+1),
     &                       Q( (iex(i+1,k)-1)*nbasis+1 ),nbasis)
6              continue
7           continue
         end if
ccc         call gmem_free_inf(kscr,'vbin.m','vbfirst','kscr')
ccc         call gmem_free_inf(ks,'vbin.m','vbfirst','ks')
      end if 

_IF(debug)
      write(iwr,*) '(_) VBFIRST - NO IF BLOCK AFTER (2) __ OK'
_ENDIF
      kvtmp = igmem_alloc_inf(kkk,'vbin.m','vbfirst','kvtmp',
     &                         IGMEM_DEBUG)

c...  nullify a temporary Vtmp vector  
      call vclr(Q(kvtmp),1,kkk)

c...  copies V vector into temporary Vtmp vector
      call dcopy(kkk,v,1,Q(kvtmp),1)

c...  nullify primary V vector
      call vclr(v,1,kkk)

      kk = 1
c
c...  move the core orbitals
c
      if (ncore.gt.0) then
_IF(debug)
      write(iwr,*) '(3) VBFIRST - NCORE > 0 __ OK'
_ENDIF
         do 10 i=1,ncore
            if (mapcie(i).gt.ncol) then
               write(iwr,*) 'i ',i,mapcie(i),' ncol ',ncol
               call vberr('core orbital index > ncol')
            end if
            call dcopy(nbasis,Q(kvtmp+(mapcie(i)-1)*nbasis),1,
     &                 v(kk),1)
            kk = kk + nbasis
10       continue
         if (logcor) then
            write(iwr,606) ncore
606         format(//' -- the ',i3,' core orbitals  --')
            call prvc(v,ncore,nbasis,v,'o','l')
         end if
      end if

_IF(debug)
      write(iwr,*) '(_) VBFIRST - NO IF BLOCK AFTER (3) __ OK'
_ENDIF
c...  move the active orbitals

      do 20 i=1,nsa+ntemp
        if (mapiee(i).gt.ncol) then
           write(iwr,*) ' ncol ',ncol,' # indices ',nsa+ntemp,
     1                  ' indices ',(mapiee(ii),ii=1,nsa+ntemp)
           call vberr('active orbital index > ncol')
        end if
        call dcopy(nbasis,Q(kvtmp+(mapiee(i)-1)*nbasis),1,
     &             v(kk),1)
        kk = kk + nbasis
20    continue
      call gmem_free_inf(kvtmp,'vbin.m','vbfirst','kvtmp')
c
c...  do mix -option on active orbitals
c
      if (nmix.gt.0) then
_IF(debug)
         write(iwr,*) '(4) VBFIRST - NMIX > 0'
_ENDIF
         nnn = (ncore+nsa+ntemp)*nbasis
         nna = ncore*nbasis
         call fmove(Q(nna+1),Q(nnn+1),(nsa+ntemp)*nbasis)
         write(iwr,608) nmix
608      format(/' -- mixing of',i3,' active orbitals requested --',/,
     *           ' orbital i . coeff +  orbital . coeff .......')
         do 50 i=1,nmix
            nnorg = nna + (imix(1,i)-1)*nbasis
            call vclr(Q(nnorg+1),1,nbasis)
            do 40 j=1,nopomx(i)
               do 30 k=1,nbasis
                  Q(nnorg+k) = Q(nnorg+k)
     *                + Q(nnn+(imix(j,i)-1)*nbasis+k) * cmix(j,i)
30             continue
40          continue
          write(iwr,609) (imix(j,i),cmix(j,i),j=1,nopomx(i))
609       format(3x,i3,f12.5,(t22,9(i3,f9.5),1x))
50        continue
      end if
c
1000  ncol = nsa + ncore
c
      kmaxmem = igmem_max_memory()
      kmemscr = kmaxmem/10
      kscr = igmem_alloc_inf(kmemscr,'vbin.m','vbfirst','kscr',
     &                       IGMEM_DEBUG)
      call bond1(Q(kscr),kmemscr)
      call gmem_free_inf(kscr,'vbin.m','vbfirst','kscr')
      if (hguess.and.ofirst) then
_IF(debug)
         write(iwr,*) '(5) VBFIRST - HGUESS and OFIRST __ OK'
_ENDIF
         kvec = igmem_alloc_inf(nbasis*nbasis,'vbin.m','vbfirst','kvec',
     &                          IGMEM_DEBUG)
         kir = igmem_alloc_inf(nbasis,'vbin.m','vbfirst','kir',
     &                         IGMEM_DEBUG)
         kic = igmem_alloc_inf(nbasis,'vbin.m','vbfirst','kic',
     &                         IGMEM_DEBUG)
         call hybmak(v(ncore*nbasis+1),Q(kvec),Q(kir),
     +               Q(kic),nbasis,nsa,qq)
         call gmem_free_inf(kic,'vbin.m','vbfirst','kic')
         call gmem_free_inf(kir,'vbin.m','vbfirst','kir')
         call gmem_free_inf(kvec,'vbin.m','vbfirst','kvec')
      end if
c
c...  adapt to harmonic basis if harmonic is required
c
      if (oharm) then
_IF(debug)
         write(iwr,*) '(6) VBFIRST - OHARM __ OK'
_ENDIF
         otran = .false. 
         write(harmvb,'(a2,i3)') 'vb',ncol
         kscr = igmem_alloc_inf(nbasis*nbasis,'vbin.m','vbfirst',
     &                          'kscr',IGMEM_DEBUG)
c
c...  get occupied orbitals in adapted basis
c
         write(iwr,'(/1x,a16,i4,a9)')'*** harmonised ',ncol,' orbitals'
         call anorm(Q(kscr),Q(1))
         call tback(v,ilifq,v,ilifq,ncol)
c...  zeroise "cartesians" to be sure
         do i=1,nbasis
            if (ielimh(i).ne.0) then
               do j=1,ncol
                  v(i+(j-1)*nbasis) = 0.0d0
               end do
            end if
         end do
c
c...  now transform back to unadapted
c
         call tdown(v,ilifq,v,ilifq,ncol)
         otran = .true.
         call gmem_free_inf(kscr,'vbin.m','vbfirst','kscr')
      end if
c
      if (logact) then
_IF(debug)
         write(iwr,*) '(7) VBFIRST - LOGACT __ OK'
_ENDIF
         write(iwr,607) ncol
607      format(//' -- the ',i3,' active+core orbitals  --')
         call  prvc(v,ncol,nbasis,v,'o','l')
      end if
c
      if (nscf.gt.0) then
_IF(debug)
         write(iwr,*) '(8) VBFIRST - NSCF __ OK'
_ENDIF
c...  finish scf preparation
c
c...  perform orthogonolisations as defined in forbex
c
         call restror(v(ncore*nbasis+1),nscf,nbasis,Q(1))

c...     clear vectors if hybrids are to be made
c
         if (hybry.and.clean) then
_IF(debug)
           write(iwr,*) '(9) VBFIRST - NSCF HYBRY CLEAN __ OK'
_ENDIF
           call clvecvb(v(ncore*nbasis+1),nscf,nbasis,'active','any',
     1                  'in vbfirst')
         end if
c
c...  next call was originally first call in subroutine equiv
c...  however, nsingly and isingly are initialised in it, and
c...  these are needed in virtual (see below) 31-08-1992
c...  this use of scr is not easy to avoid and is well protected
c
         kmaxmem = igmem_max_memory()
         kmemscr = kmaxmem/2
         kmemleft = kmaxmem - kmemscr
         kscr = igmem_alloc_inf(kmemscr,'vbin.m','vbfirst','kscr',
     &                          IGMEM_DEBUG)
         call doubly(Q(kscr),kmemscr)
         call gmem_free_inf(kscr,'vbin.m','vbfirst','kscr')
c
c...     construct the virtual space
c...
c
         ncol  = ncol + nbasis
         kvirt = 1 + (nscf+ncore)*nbasis
         maxvirt = maxvec - (nscf+ncore)
         if (maxvirt.lt.mxorbvb) write(iwr,'(a,I8)') 
     1       ' *** maxvirt in vbfirst ',maxvirt
c
c...     include core-orbitals
c
         call virtual(v,v(kvirt),nbasis,ncore,ncore+nscf,
     1                nvirt,maxvirt)
c
c...  reformat equivalence data and check consistency, and possibly
c...  automatically determine excitation patterns and equivalence
c...  relations
c
         nmo = nscf
         nao = nbasis
         naomo = nao + nmo
         nnn = naomo*(naomo+1)/2
         kaomoh = igmem_alloc_inf(2*nnn,'vbin.m','vbfirst',
     &                            'kaomoh',IGMEM_DEBUG)
         call equivt(v(ncore*nbasis+1),Q(kaomoh),nao,nmo,
     +               Q(kaomoh+nnn))
         call gmem_free_inf(kaomoh,'vbin.m','vbfirst','kaomoh')
c
         nnscf = nscf
         do i=1,nscf
            if (ofrzmo(i)) nnscf = nnscf - 1
         end do
c
         write(iwr,651) nnscf
651      format(//,
     &     13x,' -----------------------------------------------',
     *   /,13x,' =>   optimisation of',i3,' orbitals requested   <=',
     &   /,13x,' -----------------------------------------------')
c
         kk = 1
         im = 0
         if (autscf) im = nscf
         do 130 i=1,nequi
            if (iequi(i).eq.1) then
               test = '    '
            else
               write(test,'(a2,i2)') 'eq',i
            end if
            do 120 j=1,iequi(i)
               call icopy(nexw(kk)+1,iexw(1,kk),1,iex(1,kk),1)
               nex(kk) = nexw(kk)
               if (j.gt.1) test = '    '
               iexw(1,kk) = iexw(1,kk) * ieqsig(1,kk)
               do 121 k=2,nexw(kk)+1
121            iexw(k,kk) = (iexw(k,kk)-im)*iflip(iabs(iexw(1,kk)),
     &                                                 iexw(k,kk))
               if (exact(i).and.j.eq.1.and.iequi(i).ne.1) then
      if (nexw(kk).gt.1)  
     1      write(iwr,652) test,' exact',(iexw(k,kk),k=1,nexw(kk)+1)
               else if(iequi(i).ne.1.and.j.eq.1) then
      if (nexw(kk).gt.1)    
     1      write(iwr,652) test,'forced',(iexw(k,kk),k=1,nexw(kk)+1)
               else if (iequi(i).eq.1) then
      if (nexw(kk).gt.1)    
     1      write(iwr,652) test,'      ',(iexw(k,kk),k=1,nexw(kk)+1)
               else
      if (nexw(kk).gt.1)
     1      write(iwr,652) test,'     -',(iexw(k,kk),k=1,nexw(kk)+1)
               end if
c.....         iflip is equivalence sign
               iexw(1,kk) = iabs(iexw(1,kk))
               do 122 k=2,nexw(kk)+1
122            iexw(k,kk) = iabs(iexw(k,kk))+nscf
120         kk = kk + 1
130      continue
c
652      format(1x,a4,1x,a6,1x,i3,': ',(t19,25i4))
c
         if (ncurt.gt.0) write(iwr,654) ncurt
654      format(/,' 4-index curtailed to : ',3x,i4)
      end if
c
c      write(iwr,655) 
c655   format(' === The SCF orbitals ===')
c      call prvc(q(kvec),ncol,nbasis,q(kvec),'o','l')
c
c... orbitals have been rearranged change maps accordingly for
c... next point
c
      do i=1,ncore
        mapcie(i)=i
      end do
      do i=1,nsa
        mapiee(i)=i+ncore
      end do
c
c...  final check
c
_IF(atmol)
         rewind 25
         read(25) ndum1,ndum2,ndum3,ndum4,ndum5,ndum6,ndum7,
     &         ndum8,ndum9,imax
_ELSE
         call wr15vb(ndum1,ndum2,ndum3,ndum4,ndum5,ndum6,ndum7,
     &               ndum8,ndum9,imax,ndum11,ndum12,'read')
_ENDIF
         if (imax.ne.nsa) then
            write(iwr,6903) nsa,imax
6903        format(/,' The number of active orbitals is      ',i5,
     &             /,' The highest occupied orbital number is',i5,
     &             /,' This is *bad* ',/,1x)
            call vberr('# active orbitals inconsistent with crestr')
         end if
c     
      return
      end

c***********************************************************************
      subroutine atom_hybrid(v,vv,h,t,maskh,ihat,scr)
c
c...  use information from gamess and vectors just read (ncol > 0) 
c...  and the atoms specified to determine the ao's in the hybrids
c...  remember to skip core
c
      implicit REAL  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
c
      dimension v(*),vv(*),h(*),t(*),maskh(*),ihat(*),scr(*)
c       ncol*nbasis,(ncol+nbasis)*nbasis,tri(nbasis),
c       tri(nacat(natom)+nbasis),nbasis,nbasis,nbasis+nacat(natom)
c
INCLUDE(common/tractlt)
INCLUDE(../m4/common/sizes)
INCLUDE(common/turtleparam)
INCLUDE(common/infato)
INCLUDE(common/vbcri)
INCLUDE(common/splice)
INCLUDE(../m4/common/runlab)
INCLUDE(../m4/common/infoa)
c
      if (ncol.le.0) call caserr('vbvectors should precede hybrid')
c
c...  nac = # active vectors in this hybrid
c
      nac = nacat(natom)
c
      call rdistr(ihat,nhat,nbasis)
      do i=1,nbasis
         maskh(i) = 0
         if (nat.lt.100) then
            read(zbflab(i)(1:2),'(i2)') iatom
         else if (nat.ge.100.and.nat.lt.1000) then
            read(zbflab(i)(1:3),'(i3)') iatom
         else
            call caserr(' nat > 1000 not allowed in vb')
         end if
         do j=1,nhat
            if (ihat(j).eq.iatom) maskh(i) = 1
         end do
      end do
c
c...  clean mo's in this set   (i.e. clean is mandatory and might be fatal)
c
      do j=1,nac
         ibase = (iacat(j,natom)+ncore-1)*nbasis
         do i=1,nbasis
            if (maskh(i).eq.0) v(ibase+i) = 0.0d0
         end do
c...     see if we accidently deleted a complete vector
         dn = ddot(nbasis,v(ibase+1),1,v(ibase+1),1)
         if (dn.lt.1.0d-10) then
           write(iwr,*) ' *** in cleaning mo ',iacat(j,natom)+ncore,
     1                  ' for atom',iatom,' as required by atom_hybrid',
     1                  ' the mo is gone - norm ',dn,' ***'
           write(iwr,*) ' *** CHECK orbital/atom definitions ***'
           call vberr(' atom_hybrid deleted orbital ')
         end if
c
         call dcopy(nbasis,v(ibase+1),1,vv((j-1)*nbasis+1),1)
      end do
c
c...  get h-t matrix
c
      call get1e(h,dummy,'v',t)
c
c...  generate v + unit matrix and get h in this basis (=> t)
c
      call vclr(vv(nac*nbasis+1),1,nbasis*nbasis)
      do i=1,nbasis
         vv((nac+i-1)*nbasis+i) = 1.0d0
      end do
      call maket(t,vv,nac+nbasis,nbasis,h,scr)
c
c...  make a vector which indicates which ao's interact (<> 0)
c
      do i=1,nbasis
         h(i) = 0.0d0
         do j=1,nac
            h(i) = h(i) + dabs(t((nac+i)*(nac+i-1)/2+j))
         end do
         if (maskh(i).eq.0.or.h(i).lt.crihyb) h(i) = 0.0d0
      end do
c
c...  make the hybrid ao-definitions
c
      nopa(natom) = 0
      do i=1,nbasis
         if (h(i).ne.0.0d0) then
            nopa(natom) = nopa(natom) + 1
            if (nopa(natom).gt.maxopa) call caserr('maxopa overflow')
            iopa(nopa(natom),natom) = i
         end if
      end do
c
      return
      end
      subroutine atom_hybridn(v,vv,h,t,maskh,ihat,
     1                        iacat,nacat,iopa,nopa,maxnopa,scr)
c
c...  use information from gamess and vectors just read (ncol > 0) 
c...  and the atoms specified to determine the ao's in the hybrids
c...  remember to skip core
c...  general version intended to replace atom_hybrid
c...  take care to keep the ao order
c
      implicit REAL  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
c
      integer iacat(nacat),iopa(maxnopa)
      dimension v(*),vv(*),h(*),t(*),maskh(*),ihat(*),scr(*)
c       ncol*nbasis,(ncol+nbasis)*nbasis,tri(nbasis),
c       tri(nacat(natom)+nbasis),nbasis,nbasis,nbasis+nacat(natom)
c
INCLUDE(common/tractlt)
INCLUDE(../m4/common/sizes)
INCLUDE(common/turtleparam)
INCLUDE(common/vbcri)
INCLUDE(common/splice)
INCLUDE(../m4/common/runlab)
INCLUDE(../m4/common/infoa)
c
      if (ncol.le.0) call caserr('vbvectors should precede hybrid')
c
c...  nacat= # active vectors in this set / iacat numbers
c...  read atom numbers
c
      call rdistr(ihat,nhat,nbasis)
      do i=1,nbasis
         maskh(i) = 0
         if (nat.lt.100) then
            read(zbflab(i)(1:2),'(i2)') iatom
         else if (nat.ge.100.and.nat.lt.1000) then
            read(zbflab(i)(1:3),'(i3)') iatom
         else
            call caserr(' nat > 1000 not allowed in vb')
         end if
         do j=1,nhat
            if (ihat(j).eq.iatom) maskh(i) = j
         end do
      end do
c
c...  clean mo's in this set   (i.e. clean is mandatory)
c
      do j=1,nacat
         ibase = (iacat(j)+ncore-1)*nbasis
         do i=1,nbasis
            if (maskh(i).eq.0) v(ibase+i) = 0.0d0
         end do
         call dcopy(nbasis,v(ibase+1),1,vv((j-1)*nbasis+1),1)
      end do
c
c...  get h-t matrix
c
      call get1e(h,dummy,'v',t)
c
c...  generate v + unit matrix and get h in this basis (=> t)
c
      call vclr(vv(nacat*nbasis+1),1,nbasis*nbasis)
      do i=1,nbasis
         vv((nacat+i-1)*nbasis+i) = 1.0d0
      end do
      call maket(t,vv,nacat+nbasis,nbasis,h,scr)
c
c...  make a vector which indicates which ao's interact and are needed (<> 0)
c
      do i=1,nbasis
         h(i) = 0.0d0
         do j=1,nacat
            h(i) = h(i) + dabs(t((nacat+i)*(nacat+i-1)/2+j))
         end do
         if (maskh(i).eq.0.or.h(i).lt.crihyb) maskh(i) = 0
      end do
c
c...  make the hybrid ao-definitions
c
      nopa = 0
      do j=1,nhat
         do i=1,nbasis
            if (maskh(i).eq.j) then
               nopa = nopa + 1
               if (nopa.gt.maxopa) call caserr('maxopa overflow')
               iopa(nopa) = i
            end if
         end do
      end do
c
      return
      end
      subroutine deldoub(idoub,idet,nelec,nalpha,ndoub,ndets)
c
      implicit none
c
      integer nelec,nalpha,ndoub,ndets,i,j,k,jj,detno
      integer idet(nelec*ndets)
      integer idoub(ndoub)
      do i=1,ndets
         detno = (i-1)*nelec
         do j=1,nelec
            do k=1,ndoub
               if (idet(detno+j).eq.idoub(k)) idet(detno+j)=0
            enddo
            if (idet(detno+j).ne.0) then
               jj=0
               do k=1,ndoub
                  if (idoub(k).lt.idet(detno+j)) jj=jj+1
               enddo
               idet(detno+j)=idet(detno+j)-jj
            endif
         enddo
      enddo
      jj=0
      do j=1,nelec*ndets
         jj=jj+1
         if (idet(jj).eq.0) then
            do k=jj,nelec*ndets-1
               idet(k)=idet(k+1)
            enddo
            jj=jj-1
         endif
      enddo
      return
      end
