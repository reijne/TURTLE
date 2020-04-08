c
c  $Author: jvl $
c  $Date: 2014-09-15 17:36:10 +0200 (Mon, 15 Sep 2014) $
c  $Locker:  $
c  $Revision: 6297 $
c  $Source: /c/qcg/cvs/psh/GAMESS-UK/vb/vbci.m,v $
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
c  Revision 1.41  2006/01/13 17:56:47  jmht
c  Merging in changes from the release branch
c
c  Revision 1.40  2005/09/23 14:34:55  jmht
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
c  Revision 1.39.2.4  2005/11/15 11:28:55  hvd
c  Bracketting catastrophoid. Removed superflous bracket.
c
c  Revision 1.39.2.3  2005/10/28 13:43:31  jvl
c  added text to indicate from where a davidson is called
c
c  Revision 1.39.2.2  2005/07/19 06:52:19  wab
c
c  min0/max0 changes (see m4 checkin)
c
c  Revision 1.39.2.1  2005/06/24 05:42:34  wab
c  Second attempt to do this after the bug in guess that caused apparent
c  failures to be diagnosed yesterday has been fixed
c
c  Revision 1.39  2005/04/22 11:07:52  jvl
c  All vb criteria are now in common vbcri and may be set
c  in the input. Also their defaults are now consistent
c
c  Revision 1.38  2005/03/31 14:59:58  jvl
c  - corrected pert option, which was screwed earlier
c  - include file hsinfo
c
c  Revision 1.37  2005/03/25 23:31:32  jvl
c  changes to allow i8 VB; labeling of CI matrix elements => packed real
c
c  Revision 1.36  2005/02/05 18:04:24  wab
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
c  1. all explicit specifications of "real*8" replaced with real*8 in line
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
c  Revision 1.35  2005/01/05 10:22:52  jvl
c  Bugfix concerning erasing scratch information between consequetive VBSCF
c  runs (for PCM). Further extension of the PCM model with a tweak input
c  option (JJE).
c
c  Revision 1.34  2004/07/15 15:55:09  jvl
c  Dimensioning bug fixed (JJE)
c
c  Revision 1.33  2004/05/19 22:57:44  jvl
c  Brought harmonic to VB; If hamonic is specified the cartesion
c  components are projected out of the occupied and virtual orbitals
c  yielding a harmonic vbscf and vbci (the latter may of may not be what you want)
c
c  Revision 1.32  2004/04/26 11:05:05  jvl
c  Bugfix concerning printout of more than 240 weights (JJE)
c
c  Revision 1.31  2003/11/19 16:50:23  jvl
c  bugfix concerning confusion with coeff and bcoeff, resp. coefficients for
c  determiniants in structures and coefficients of determinants in psi0 (JJE).
c
c  Revision 1.30  2003/08/26 13:35:56  jvl
c  Weights are printed together with their branching or rumer diagram.
c  (JJE 26/8/2003)
c
c  Revision 1.29  2003/08/20 09:49:36  jvl
c  bugfix regarding wrong dimension of igroup array (JJE 20/8/2003)
c
c  Revision 1.28  2003/06/10 16:40:41  jvl
c  davdia was expected to converge to the eigenvalue with the biggest coefficient
c  on position 1 (psi0 in bi-vector). This was not the true in the case of li2.
c  Changed it (JJE 10/6/2003).
c
c  Revision 1.27  2003/03/21 14:08:23  jvl
c  Fixed a bug concerning dimensioning an array (iorth). JJE 21-3-2003
c
c  Revision 1.26  2003/02/28 10:58:32  jvl
c  Fixed bug in subroutine makepsi0. Lowered vbscf.m deck from -O3 to -O2
c  for sgi systems (JJE).
c
c  Revision 1.25  2003/02/18 17:17:12  jvl
c  A lot of changes, in general:
c  - In vbci the call to bsmat for natural orbitals is no longer necessary
c  - The suboptions pert and fock for optimise in scf are implemented
c  - The subroutine vbscf is completely restructured and is more logical
c  - Addition of the makepsi0, which uses a small extra 2-electron,
c    transformation
c    JJE - 18 feb 2003
c
c  Revision 1.24  2003/01/19 12:33:06  jvl
c  - check on # active orbitals (an infinite loop could result)
c  - hcpu option to allow primitive restarting of matrix element calculation in vbci
c
c  Revision 1.23  2002/11/25 15:25:30  jvl
c  There was a restriction to 500000 words for 4-index. Not very relevant
c  and causing multipass mode, which unfortunately crashes in parallel runs.
c  Eliminated restriction and added error-message for multi-pass in parallel
c
c  Revision 1.22  2002/09/05 14:20:45  jvl
c  Print of perturbation h- and s-matrix (2 columns each) added to the
c  normal print of h- and s-matrix per brilliounstructure (JJE)
c
c  Revision 1.21  2002/05/28 15:07:46  jvl
c  General cleanup
c  got rid of all fortran-files (now part of ed7 is used)
c  handle 3 more common's through include
c  make vector sections in vb more logical
c  added print parameter to getqvb
c  ....
c
c  Revision 1.20  2002/05/16 16:20:34  jvl
c  davidson was forced (by me) to do a cycle too many ....
c
c  Revision 1.19  2002/04/28 22:30:28  jvl
c  Added virtual options for vb-orbital optimisation
c  in the process an adapted version of poporb is made
c
c  Revision 1.18  2002/04/25 17:21:33  jvl
c  added checks on negative norm on vector in davidson (= overcomplete sywstem)
c  changed a few (print)flags
c
c  Revision 1.17  2002/02/10 23:53:28  jvl
c  forgot to protect the vb-vector in call to prgsmat  ; fixed
c
c  Revision 1.16  2002/02/10 21:06:59  jvl
c  cleaned up printing and made comparison ci/scf possible
c
c  Revision 1.15  2002/02/07 12:33:31  jvl
c  fixed bugs + added print of orthog h-matrix
c
c  Revision 1.14  2002/01/17 08:49:34  jvl
c  Fixed problems with zero norm states in parallel version;
c  Was more difficult to spot in calculation stage, since some processors
c  may not calculate the diagonal amtrix element at all, zo 0.0 does not mean that much
c  then.
c
c  Revision 1.13  2001/11/28 15:18:41  jvl
c  bug in vbci; it wrote the weights to unit 7 instead of iwr; fixed
c
c  Revision 1.12  2001/10/15 15:32:39  jvl
c  major parallel overhaul; makes asynchronous switchable
c  no luck yet however ...........
c
c  Revision 1.11  2001/09/19 16:02:42  jvl
c  some corrections resulting from DEC compile
c
c  Revision 1.10  2001/07/03 14:13:15  jvl
c  fixed a few bugs in parallel code and allowed non-usage of peigss in vb
c  (i.e. peigss specifies if peigas is to be used in vb)
c
c  Revision 1.9  2001/06/28 09:06:45  jvl
c  Jacobi davidson checkin voor VB (generalised eigenvalue problem)
c  Xiny 2001
c
c  Revision 1.8  2001/06/27 15:28:46  jvl
c  Changed vector printing
c  fixed various bugs (related to getin2)
c  added frozen en frozen hybrids options (experimental)
c
c  Revision 1.7  2001/06/21 16:38:19  jvl
c  added super hybrid option (experimental and perhaps useless) and a bit of cleaning up
c
c  Revision 1.6  2001/06/15 15:16:16  jvl
c  added extra in veccomb (vb) + cleanup
c
c  Revision 1.5  2001/06/12 12:21:16  jvl
c  - Removed several bugs
c  - Improved parallel 4-index transformation including changes in getin2
c  - Added timing for VB routines
c
c  Revision 1.4  2001/06/04 21:59:19  jvl
c  added vectors colmbine option + cleaned vb print
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
c  Revision 1.1  2000/02/16 12:20:41  jvl
c  adding vb files to gamess repository
c
c Revision 1.11  1998/02/16  14:22:44  fokke
c added counter process to parallel implementation
c added output of energies and interactions of single determinants
c
c Revision 1.11  1998/02/16  14:22:44  fokke
c added counter process to parallel implementation
c added output of energies and interactions of single determinants
c
c Revision 1.10  1997/05/30  11:10:21  fokke
c Changed icopy to icopy and removed variables icopy
c
c Revision 1.9  1997/05/22  12:53:09  joop
c changed imove to icopy
c
c Revision 1.9  1997/05/22  12:53:09  joop
c changed imove to icopy
c
c Revision 1.8  1997/05/22  11:31:48  joop
c Added include macro.cpp
c
c Revision 1.7  1997/05/21  11:36:03  joop
c some fixes (21-5-97) FD
c
c Revision 1.6  1997/01/02  16:21:22  joop
c MPI changes for SGI + 1 bugfix
c
c Revision 1.5  1996/11/07  17:09:15  joop
c changed checking on degeneracies in FEW and adapted core partitioning
c
c Revision 1.4  1996/10/29  14:12:57  joop
c parallel/labeled version
c
c Revision 1.3  1996/10/29  13:28:34  joop
c added rcs comments
c
c
      subroutine davdia(length,nvec,dh,ds,cc,yc,zc,c,y,z,
     *                  vnrm,hsmall,hsave,varia,q,scr,eigen,potnuc,
     *                  energy,caller)
      implicit real*8  (a-h,o-z) ,  integer   (i-n)
c
      dimension dh(length),ds(length),cc(length),yc(length),zc(length),
     *          c(length,nvec),y(length,nvec),z(length,nvec),
     *          vnrm(*),hsmall(*),hsave(*),varia(3,*),q(nvec,nvec),
     *          scr(*),eigen(*)
      character*(*) caller
c.....
c.....non-orthogonal davidson diagonaliser
c..... with emin lock and vmin modes
c.....
c..... dh (diagonal of H) is already normalised upon entry
c..... ds (diagonal of S) is changed to 1.0d0/dsqrt(ds) here
c..... vmin needs auxilary array varia(3,nvec*nvec+1/2) :
c....     varia : (1,..) = zz; (2,..) = yz+zy; (3,..) = yy
c....  Joop van Lenthe // juli 1990
c.....
      logical ohmat
      integer nstruc,iprinv,ihfile,ihbloc,ipbloc,idumhs
      common/hsinfo/ohmat(3),nstruc,iprinv,ihfile,ihbloc,ipbloc,idumhs
      common/davcom/eshift,thresj,thresd,maxcyc,maxdav,maxsel,iprind,
     *              alter,firsts
      common/posit/iky(3)
      common/davctx/mode
c
      
      integer ird, iwr, ipu, main, ibl2m, idaf, ibl3d, num8
      integer lenbl, ntapes, iblkci, numlib, iblkpl
      integer numdis, ibldis, irdcpy
      common /iofile/ ird,iwr,ipu,main,ibl2m,idaf,ibl3d,num8,
     +                lenbl(4),ntapes(5),iblkci,numlib,iblkpl,
     +                numdis,ibldis,irdcpy
c
      character*4 mode
      logical alter,firsts,print
      logical first1
      first1 = .false.
      if (dabs(dabs(cc(1))-1.0d0).lt.thresj) first1 = .true.
c
c...  adjust ds
c
      do 1 i=1,nstruc
1     ds(i) = 1.0d0/dsqrt(ds(i))
      print = iprind.ge.10
      tester = 0.0d0
c.....
c.....form hc ( = z ) and sc ( = y )
c.....
      call hsc(zc,yc,cc,nstruc,ihfile,ihbloc,'b')
      call fmove(cc,c(1,1),nstruc)
      call fmove(yc,y(1,1),nstruc)
      call fmove(zc,z(1,1),nstruc)
c.....
c.....determine the trial energy (chc/csc) and save it
c.....
      ncycle = 1
9999  continue 
      if (ncycle.gt.1) write(iwr,99) ncycle
99    format(' '/,' davidson folded at cycle ',i4)
      if (print) write(iwr,88)
88    format(' '/,'   cycle                 energy            tester',
     &       '           variance')
      vnrm(1) = ddot(nstruc,c,1,y,1)
      hsave(1) = ddot(nstruc,c,1,z,1)
      vnrm(1) = 1.0d0 / vnrm(1)
      hsave(1) = hsave(1) * vnrm(1)
      vnrm(1) = dsqrt(vnrm(1))
      h11      = hsave(1)
c.....
c.....calculate the variance-auxilaries(if needed), and save it
c.....
      if (mode.eq.'vmin') then
          varia(1,1) = ddot(nstruc,z,1,z,1)
          varia(2,1) = 2.0d0*ddot(nstruc,y,1,z,1)
          varia(3,1) = ddot(nstruc,y,1,y,1)
      end if 
c.....
c.....start looping over nvec
c.....
      do 170 i=2,nvec
         if (maxcyc-ncycle.eq.10.and..not.print) write(iwr,88)
         print = iprind.ge.10.or.(maxcyc-ncycle).le.10  
c
         var = ddot(nstruc,zc,1,zc,1)
     *         - 2 * h11 * ddot(nstruc,zc,1,yc,1) +
     &         h11 * h11 * ddot(nstruc,yc,1,yc,1)
c.....
	 tester = 0.0d0
         do 10 j=1,nstruc
            c(j,i) = ( h11  *  yc(j)  - zc(j) ) * ds(j)
            diago  = dh(j) - (h11-eshift)
            if (dabs(diago).le.0.05d0) diago = dsign(0.05d0,diago)
            c(j,i) = c(j,i) / diago
            if (dabs(c(j,i)).ge.0.2d0) c(j,i) = dsign(0.2d0,c(j,i))
            if (dabs(c(j,i)).gt.tester) tester = dabs(c(j,i))
            c(j,i) = c(j,i)*ds(j)
10       continue
         ncycle = ncycle + 1
         if (print) write(iwr,11) ncycle,h11+potnuc,tester,var
11       format(i8,e23.14,e18.9,e19.10)
         if (h11.gt.dh(1)*(1-1d-14)) then
            write(iwr,'(a,2f25.15)') 'davidson h11 dh ',h11,dh(1) 
            call caserr('Davidson gone mad')
         end if
c...      we want to establish convergence but don't go outside space
         if (tester.lt.thresd.and.ncycle.ge.2.or.ncycle.gt.nstruc) then
            goto 190
         end if
         if (ncycle.gt.maxcyc) then
            goto 200
         end if
         if (alter) eshift = -eshift
c.....
c.....   project current space out of new c-vector (orthogonalise)
c.....
         call hsc(dum,y(1,i),c(1,i),nstruc,ihfile,ihbloc,'s')
         x  = ddot(nstruc,y(1,i),1,c(1,i),1)
         iortim = 0
20       xx = x
         if (x.lt.0.0d0) then
            write(iwr,21) x,i,tester,caller
21          format(' Negative norm detected : ',e12.5,' at cycle',i3,
     1             ' tester ',e12.5,' finishing davidson, good luck',/,
     2             ' You might consider the super option',/,
     3             ' ** used for    ',a10,' **') 
            go to 200
         end if
         iortim = iortim + 1
         do 40 j=1,i-1
            x = ddot(nstruc,c(1,i),1,y(1,j),1) * vnrm(j) * vnrm(j)
            do 30 m=1,nstruc
               c(m,i) = c(m,i) - x * c(m,j)
30          continue
40       continue
         call hsc(dum,y(1,i),c(1,i),nstruc,ihfile,ihbloc,'s')
         x = ddot(nstruc,c(1,i),1,y(1,i),1)
         if (x.lt.0.05d0*xx.and.iortim.lt.10) then
c.....
c.....      orthogonalise again
c.....
            goto 20
         end if
c
         x = 1.0 d0/ x
         vnrm(i) = dsqrt(x)
         call hsc(z(1,i),dum,c(1,i),nstruc,ihfile,ihbloc,'h')
c
         do 50 j=1,i
            temp = vnrm(i) * vnrm(j)
            hsave(iky(i)+j) = temp * ddot(nstruc,c(1,i),1,z(1,j),1)
            if (mode.eq.'vmin') then
               varia(1,iky(i)+j) = temp * ddot(nstruc,z(1,i),1,z(1,j),1)
               varia(2,iky(i)+j) = temp * 
     &   (ddot(nstruc,y(1,i),1,z(1,j),1)+ddot(nstruc,z(1,i),1,y(1,j),1))
               varia(3,iky(i)+j) = temp * ddot(nstruc,y(1,i),1,y(1,j),1)
            end if
50       continue
         if (mode.eq.'lock') then
            call fmove(hsave,hsmall,iky(i)+i)
            call jacobt(hsmall,iky,i,q,nvec,eigen,2,1,thresj,scr)
            if (first1) then
c...        make sure first element is largest
               iimax = 1
               qqmax = dabs(q(1,1))
               do k=2,i
                  if (dabs(q(1,k)).gt. qqmax) then
                     iimax = k
                     qqmax = dabs(q(1,k))
                  end if
               end do
               if (iimax.ne.1) then
                  write(iwr,'(a,i4,a)') 'Warning: Taking vec',iimax,
     &                     ' instead of 1 in davdia'
                  write(iwr,'(a,9f12.5)') ' e0+eigenvalues :',hsave(1),
     &                                    (eigen(l),l=1,min(8,iimax))
                  write(iwr,'(a,8f12.5)') ' vec(1,..) :',
     &                                    (q(1,l),l=1,min(8,iimax))
                  do k=1,i
                     q(k,1) = q(k,iimax)
                  end do
                  eigen(1) = eigen(iimax)
               end if
            end if
            if (iprind.gt.15) write(iwr,52) eigen(1) + potnuc
52          format(/' ','first eigenvalue of jacobi subspace problem :'
     &            ,e23.14)
         else if (mode.eq.'emin') then
            call fmove(hsave,hsmall,iky(i)+i)
            call jacobt(hsmall,iky,i,q,nvec,eigen,2,2,thresj,scr)
            if (iprind.gt.15) write(iwr,53) eigen(1) + potnuc
53          format(/' ','lowest eigenvalue of jacobi subspace problem :'
     &            ,e23.14)
         else if (mode.eq.'vmin') then
            do 70 j=1,iky(i+1)
               hsmall(j) = varia(1,j)-h11*varia(2,j)+h11*h11*varia(3,j)
70          continue
            call jacobt(hsmall,iky,i,q,nvec,eigen,2,1,thresj,scr)
         end if
c
         call vclr(cc,1,nstruc)
         call vclr(zc,1,nstruc)
         call vclr(yc,1,nstruc)
         do 130 m=1,i
            cexpan = q(m,1) * vnrm(m)
            call daxpy(nstruc,cexpan,c(1,m),1,cc,1)
            call daxpy(nstruc,cexpan,z(1,m),1,zc,1)
            call daxpy(nstruc,cexpan,y(1,m),1,yc,1)
120         continue
130      continue
         h11 = ddot(nstruc,cc,1,zc,1) / ddot(nstruc,cc,1,yc,1)
170   continue
c
c...   max. number of expansion-vectors reached  .. fold
c...   recalculate everything    
c
      call hsc(dum,yc,cc,nstruc,ihfile,ihbloc,'s')
      cnorm = 1.0d0/dsqrt( ddot(nstruc,cc,1,yc,1) ) 
      call dscal(nstruc,cnorm,cc,1) 
      call hsc(zc,yc,cc,nstruc,ihfile,ihbloc,'b')
      call fmove(cc,c,nstruc)
      call fmove(yc,y,nstruc)
      call fmove(zc,z,nstruc)
      goto 9999
c
c...   convergence
c
190   energy = h11+potnuc
      if (iprind.ge.10) write(iwr,22) ncycle
22    format(//20x,26('*')//20x,'convergence at cycle',i6//20x,26('*'))
c
      return
c
200   write(iwr,33)
33    format(//20x,7('*')//20x,'no convergence'//20x,7('*'))
      energy = h11+potnuc
c
      return
      end
      subroutine davidt(q,lword,potnuc,energy,caller)
c
      implicit real*8  (a-h,o-z) , integer   (i-n)
c
c
c...   calling and core-partitioning for davidson
c...   see that ci-vector is first in q
c
      dimension q(lword)
      character*(*) caller
c
      common/davcom/eshift,thresj,thresd,maxcyc,maxdav,maxsel,iprind,
     *              alter,firsts
      logical alter,firsts
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
c...   turtle parameters
c
      integer melec,ncorq,maxex,maxact,mbonds,maxato,
     &        maxopa,ncases,maxcon,mmsocc,mxorbvb,maxspvb
c
      parameter (melec=300)
      parameter (ncorq=2000000)
      parameter (maxex=20000,maxact=150,mbonds=150,
     &           maxato=750,maxopa=1500)
      parameter (ncases=35)
c      ... maxcon in vbcrestr ...
c92      parameter (maxcon = 100)
      parameter (maxcon = 3000)
c      ... msocc max. # singly occuppieds ; was always enough (2009) ...
      parameter (mmsocc = 30)
c..   maxspvb  maximum # spinfunctions per configuration 
c..   (42 for 10 singlys and the ultimate answer)
      parameter (maxspvb=42) 
c
c     mxorbvb may be as large as allowed by n_8_16 (32767)
c
      parameter (mxorbvb=1000)
      common/vbpert/crivarpert(4),shiftp,fpdoc,fpuoc,
     &              nvarp,ifirsp,idelp,npert,
     &              ntype,itype(mxorbvb),
     &              minfock,maxfock,minpert,maxpert,minvar,maxvarp,
     &              minfockd,maxfockd,
     &              fckp,ovbper,fockdiagonal
      INTEGER nvarp,ifirsp,idelp,npert,ntype,itype
      INTEGER minfock,maxfock,minpert,maxpert,minvar,maxvarp
      INTEGER minfockd,maxfockd
      logical ovbper,fckp,fockdiagonal
      real*8 crivarpert,shiftp,fpdoc,fpuoc
c description:
c   ntype -> 1 - specify optimisation per brillioun structure
c            2 - specify optimisation per orbital
c            3 - specify optimisation per kind (type excitation)
c   itype ->
c
c          ntype = 1
c            itype(1..maxex) - optimisation type for each (brillouin) structure
c
c          ntype = 2
c            itype(1..mxorbvb) - optimisation type for each orbital
c            
c          ntype = 3 
c            itype(1) - optimisation type for doubly to variably occupied
c            itype(2) - optimisation type for doubly to unoccupied
c            itype(3) - optimisation type for variably to variably occupied
c            itype(4) - optimisation type for variably to unoccupied
c            
c          ntype = 4
c            itype(1..mxorbvb) - optimisation type for each orbital (to be
c                               specified by clsbril
c
c   contents of itype array:
c            itype(x)=1 for Super CI optimalisation
c            itype(x)=2 for perturbation theory optimalisation
c            itype(x)=4 perturbation theory  determines choice between Super and Pert
c
c   crivarpert ~: critaria for itype=4 case per class if applicable

      common /brill/ nex(maxact) ,iex (maxex,maxact),
     &               nexw(maxact),iexw(maxex,maxact),
     &                            ieqsig(maxex,maxact),
     &               iequi(maxact),nequi,iredund(mxorbvb),nredund,
     &               inteq(maxact+1,maxact),ninteq(maxact),ninter,
     &               nscf,nullity,ortscf,schmi,reduce,autscf,
     &               exact(maxact),nbrill,idequi,
     &               ieqs(maxact),nexeq(maxact),iexeq(maxact*maxact)
      logical ortscf,schmi,reduce,autscf,exact
      integer nequi,nredund,ninter,nscf,nullity,nbrill,nex,nexw,
     &        iex,iexw,ieqsig,iequi,iredund,inteq,ninteq,idequi,
     &        ieqs,nexeq,iexeq
      logical ohmat
      integer nstruc,iprinv,ihfile,ihbloc,ipbloc,idumhs
      common/hsinfo/ohmat(3),nstruc,iprinv,ihfile,ihbloc,ipbloc,idumhs
c
c     core :  frozen core energy for the whole job
c     core_i :  frozen core energy  in makeps0 
c
      real*8 core,pot_nuc,core_i
      integer nfil,nofil,iblvb,lblvb,iblhs
      common /ffile/ core,core_i,pot_nuc,nfil,
     1               nofil(20),iblvb(20),lblvb(20),iblhs
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
c...            symcri = critor*100 remcri = 1.0d-10 crihyb = 1.0d-14
c...            cribri = 1.0d-8 crinor = critor*100  cribub = 1.0d-10
c...            criign = 1.0d-44 cripop 1.0d-10
      real*8 critor,crilow,cridep,symcri,remcri,crihyb,cribri,
     1     crinor,cribub,criign,cripop
      common/vbcri/ critor,crilow,cridep,symcri,remcri,crihyb,cribri,
     1              crinor,cribub,criign,cripop
      common/discc/ ied(16)
      character*4 ied
      logical oroot
c
      
      integer ird, iwr, ipu, main, ibl2m, idaf, ibl3d, num8
      integer lenbl, ntapes, iblkci, numlib, iblkpl
      integer numdis, ibldis, irdcpy
      common /iofile/ ird,iwr,ipu,main,ibl2m,idaf,ibl3d,num8,
     +                lenbl(4),ntapes(5),iblkci,numlib,iblkpl,
     +                numdis,ibldis,irdcpy
c
c
c...   print h and s-matrices and other info if requested via iprind
c
      if (iprind.ge.100) write(iwr,601) nstruc,core,
     *                 ied(ihfile),ihbloc,caller
601   format(/' -- enter davidson with dimension ',i10,/,
     *        ' frozen core energy ',f24.16,
     *        '    h/s-matrix on ',a4,' at block ',i7,/,
     #        ' used for ',a10)
       if ((iprind.ge.100).or.ohmat(3)) then
         nstrucc = dsqrt(lword*1.0d0)
         nstrucc = min(nstrucc,nstruc)
         write(iwr,602) 'h-matrix'
         call prhs(q,nstrucc,ihfile,ihbloc,'h')
         write(iwr,602) 's-matrix'
         call prhs(q,nstrucc,ihfile,ihbloc,'s')
         kscr=1
         if (ipbloc.gt.0) then
            write(iwr,602) 'pert-h'
            call rdedx(q(kscr),nbrill,ipbloc,ihfile)
            call reads(q(kscr+nbrill),nbrill,ihfile)
            do i=0,nbrill-1
               write(iwr,603) q(kscr+i),q(kscr+nbrill+i)
            enddo
            call reads(q(kscr),nbrill,ihfile)
            call reads(q(kscr+nbrill),nbrill,ihfile)
            write(iwr,602) 'pert-s'
            do i=0,nbrill-1
               write(iwr,603) q(kscr+i),q(kscr+nbrill+i)
            enddo
         endif
602   format(//,'  ----  the ',a8,'  ----',/)
603   format('  ',2f22.16)
       end if
c
c...   startup
c
      maxsel = min(maxsel,nstruc)
      max100 = max(maxsel,100)
      kcc = 1
      kdh = kcc + nstruc
      kds = kdh + nstruc
      ke = kds + nstruc
      kimap = ke + max100
      kimapc = kimap + max100*2
      khs = kimapc + max100*2
      khso = khs + max100*(max100+1)/2
      kq = khso + max100*(max100+1)/2
      kqt = kq + max100*max100
      kx = kqt + max100*max100
      kitrou = kx + 2*max100*max100
      kend = kitrou + nstruc
      if (kend.gt.lword) call corfait(kend,lword,'david before few')
c
c...  get h and s diagonals
c
      call getd(q(kdh),q(kds),nstruc,ihfile,ihbloc)
c
c...  eliminate zero states if not possible before (e.g.parallel)
c...  *before* normalising h-diagonal
c
      do 20 i=1,nstruc
         if (q(kds+i-1).lt.crinor) then
           if (oroot()) then
            write(iwr,2) i,q(kds+i-1)
2           format(/'              ***warning***',/
     +              ' diagonal s-element ',i9,' too small :',e12.5,/,
     +              ' element set to 1.0d0 => state removed ',/,
     +              ' run serial version to get more info ',/,
     +              ' ************************************************')
           end if
           q(kds+i-1) = 1.0d0
         end if
c...   normalise h-diagonal
         q(kdh+i-1) = q(kdh+i-1) / q(kds+i-1)
20    continue
c
      call few(q(kcc),q(kdh),q(kds),maxsel,
     *         q(ke),q(kimap),q(khs),q(khso),q(kq),q(kqt),q(kx),potnuc,
     *                                     q(kimapc),q(kitrou))
c
c...   determine maximum  expansion-vectors
c
      nn = 6*maxdav**2  + 5*maxdav + 5*nstruc
      nword = lword - nn
      maxv = nword/(nstruc*3)
      maxv = min(maxv,maxdav)
      if (maxv.le.2) call vberr(' insufficient room for davidson')
c
c...  davidson core-partioning
c
      kyc = kds + nstruc
      kzc = kyc + nstruc
      kc = kzc + nstruc
      ky = kc + nstruc*maxv
      kz = ky + nstruc*maxv
      kvnorm = kz + nstruc*maxv
      khsmall = kvnorm + maxv
      khsave = khsmall + maxv*(maxv+1)/2
      kvaria = khsave + maxv*(maxv+1)/2
      kq     = kvaria + 3 * maxv*(maxv+1)/2
      kscr   = kq + maxv*maxv
      keigen = kscr + 2*maxv
      kend = keigen + maxv
      if (kend.gt.lword) call corfait(kend,lword,' stupid david')
c
      timed = cpulft(1)
      call davdia(nstruc,maxv,
     *            q(kdh),q(kds),q(kcc),q(kyc),q(kzc),q(kc),q(ky),q(kz),
     *            q(kvnorm),q(khsmall),q(khsave),q(kvaria),q(kq),
     *            q(kscr),q(keigen),potnuc,energy,caller)
      timed = cpulft(1) - timed
      if (iprind.gt.100) write(iwr,77) timed
77    format(/' ','the diagonalisation costs ',f8.3,' seconds')
c
      return
      end

      subroutine davids(q,lword,potnuc,energy,caller)
c
c.......................................................................
c.......................................................................
c     begin of diagonalisation part
c.......................................................................
c.......................................................................
c
      implicit real*8  (a-h,o-z) ,  integer   (i-n)
c
      dimension q(lword)
      character*(*) caller
c
c...  davidson for bi-diagonalisation
c...  switches from davcom to davscf
c
      common/davcom/dav(3),idav(4),ldav(2)
      common/davctx/mode
c
      common/davscf/davs(3),idavs(4),ldavs(2)
      common/davssx/modess
c
      dimension sav(3),isav(4),lsav(2)
c
      logical ldav,ldavs,lsav
      character*4 mode,modess,smode
c
c...  switch parameters
c
      do 1  i=1,3
         sav(i) = dav(i)
         dav(i) = davs(i)
1     continue
      do 2  i=1,4
         isav(i) = idav(i)
         idav(i) = idavs(i)
2     continue
      do 3 i=1,2
         lsav(i) = ldav(i)
         ldav(i) = ldavs(i)
3     continue
      smode = mode
      mode = modess
c
c...   now call davidson or jacobi davidson
c
      if (mode.eq.'jacd') then
c
         call cjacdav(q,lword,potnuc,energy)
      else
c
         call davidt(q,lword,potnuc,energy,caller)
      end if
c
c...   restore common  /davcom/
c
      do 11  i=1,3
11    dav(i) = sav(i)
      do 22 i=1,4
22    idav(i) = isav(i)
      do 33 i=1,2
33    ldav(i) = lsav(i)
      mode = smode
c
      return
      end
      subroutine few(c,dh,ds,nsel,
     &               e,imap,hs,hso,q,qt,x,potnuc,imapc,itroub)
c
      implicit real*8 (a-h,o-z), integer (i-n)
c
c...  this routine makes up a start for the non-orthogonal davidson
c
c...  q,qt :   nsel**2
c...  x    : 2*nsel**2
c...  imap : 2*nsel
c
c...  itroub (returned by schmidt) indicates that the subspace contains
c...  dependencies. if so, the procedure is repeated without the redun-
c...  dancy. in the next iteration the corresponding excitation will be
c...  omitted
c
c
c...   turtle parameters
c
      integer melec,ncorq,maxex,maxact,mbonds,maxato,
     &        maxopa,ncases,maxcon,mmsocc,mxorbvb,maxspvb
c
      parameter (melec=300)
      parameter (ncorq=2000000)
      parameter (maxex=20000,maxact=150,mbonds=150,
     &           maxato=750,maxopa=1500)
      parameter (ncases=35)
c      ... maxcon in vbcrestr ...
c92      parameter (maxcon = 100)
      parameter (maxcon = 3000)
c      ... msocc max. # singly occuppieds ; was always enough (2009) ...
      parameter (mmsocc = 30)
c..   maxspvb  maximum # spinfunctions per configuration 
c..   (42 for 10 singlys and the ultimate answer)
      parameter (maxspvb=42) 
c
c     mxorbvb may be as large as allowed by n_8_16 (32767)
c
      parameter (mxorbvb=1000)
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
c     twice contains info on pauli-principle related restrictions
c     ("doubly","redorb","iflip"), and super, which indicates wether
c     internal excitations are allowed for (not so if super-mcscf)
c     alos ifreez/nfreez orbitals to be frozen
      integer inforb,ninfor,idoubly,ndoubly,isingly,nsingly,kaart,
     &        ioror,igno_sh,igno_sel
      logical super,super_hybrid,ofrzmo,super_cas,super_act
      common /twice/ inforb(maxact*5,maxact),ninfor(maxact),
     &               idoubly(maxact),ndoubly,isingly(maxact),nsingly,
     &               kaart(mxorbvb+maxact),ioror(maxact,maxact),super,
     &               super_hybrid,igno_sh(maxato),ofrzmo(maxact),
     &               super_cas,super_act,igno_sel
c...
      common /davcom/ eshift,thresj,thresd,maxcyc,maxdav,maxsel,iprind,
     &                alter,firsts
      logical alter,firsts
      common /forbex/ nforbi,iforbi(maxact*maxact,2),
     &                nionj,ionj(maxact,2),
     &                nsetor,isetor(maxact+1,maxact),
     &                nionres,ionrest(maxact)
      common /posit/ iky(3)
      logical ohmat
      integer nstruc,iprinv,ihfile,ihbloc,ipbloc,idumhs
      common/hsinfo/ohmat(3),nstruc,iprinv,ihfile,ihbloc,ipbloc,idumhs
      common /brill/ nex(maxact) ,iex (maxex,maxact),
     &               nexw(maxact),iexw(maxex,maxact),
     &                            ieqsig(maxex,maxact),
     &               iequi(maxact),nequi,iredund(mxorbvb),nredund,
     &               inteq(maxact+1,maxact),ninteq(maxact),ninter,
     &               nscf,nullity,ortscf,schmi,reduce,autscf,
     &               exact(maxact),nbrill,idequi,
     &               ieqs(maxact),nexeq(maxact),iexeq(maxact*maxact)
      logical ortscf,schmi,reduce,autscf,exact
      integer nequi,nredund,ninter,nscf,nullity,nbrill,nex,nexw,
     &        iex,iexw,ieqsig,iequi,iredund,inteq,ninteq,idequi,
     &        ieqs,nexeq,iexeq
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
c...            symcri = critor*100 remcri = 1.0d-10 crihyb = 1.0d-14
c...            cribri = 1.0d-8 crinor = critor*100  cribub = 1.0d-10
c...            criign = 1.0d-44 cripop 1.0d-10
      real*8 critor,crilow,cridep,symcri,remcri,crihyb,cribri,
     1     crinor,cribub,criign,cripop
      common/vbcri/ critor,crilow,cridep,symcri,remcri,crihyb,cribri,
     1              crinor,cribub,criign,cripop
c
      dimension c(*),dh(*),ds(*),e(*),imap(*),hs(*),hso(*),q(*),qt(*)
     &         ,x(*),itroub(*) ,imapc(*)
c
      logical reduct
c
      common /txema/ struc,ivec
      logical struc
c
1     nsel = min(nsel,nstruc)
      do 5 i=1,nsel
5     imap(i) = i
      max = nsel
      dhmax = dh(max)
      if (.not.firsts) then
         kk = nsel
c
10       max = 1
         do 20 i=1,nsel
            if (dh(imap(i)).gt.dh(imap(max))) max = i
20       continue
c
c...     check rest of diagonal
c
         dhmax = dh(imap(max))
         do 40 i=kk+1,nstruc
            if (dh(i).lt.dhmax) then
               imap(max) = i
               kk = i
               go to 10
            end if
40       continue
      end if
c
c...  check if there is a degeneracy
c
      reduct = .false.
      if (nsel.eq.nstruc) go to 60
      do 50 i=1,nstruc
         if (dabs(dh(i)-dhmax).lt.cribub) then
            if (i.ne.max) then
c...   check if i in set imap
               do 45 j=1,nsel
45             if (i.eq.imap(j)) go to 50
c...   sorry degenary outside choosen set => reduce set
               reduct = .true.
               go to 60
            end if
         end if
50    continue
60    nn = nsel
      if (reduct) then
c...   reduce choosen set until no degeneracy with outside world
         kk = nsel
         do 70 i=1,nsel
70       imap(kk+i) = imap(i)
         nn = 0
         do 80 i=1,nsel
            if (dabs(dh(imap(kk+i))-dhmax).gt.cribub) then
               nn = nn + 1
               imap(nn) = imap(kk+i)
            end if
80       continue
      end if
c
      if (nn.le.0) then
c...   we did eliminate all => make nsel bigger
         nsel = nsel+1
         if (nsel.gt.100) call vberr(' over 100 degenerate states ')
         go to 1
      end if
c
c...    now diagonalise choosen subspace
c
      call ibubbl(imap,nn,1)
c
c...    now get s submatrix    and determine orthogonalising vectors
c
      ntroub = 0
12345 call submat(hs,nn,imap,ihfile,ihbloc,nstruc,'s')
c      write (iwr,'(A)') ' s-submatrix'
c      call tripri(hs,nn)
c.....save submatrix for print concerning dependencies
      call fmove(hs,hso,nn*(nn+1)/2)
c.... 
      call schmidt(hs,q,nn,ktroub)
      if (ktroub.ne.0) then
c....
c....    remove the trouble-maker and try again
         ntroub = ntroub + 1
         itroub(ntroub) = imap(ktroub)
         if (nn.gt.ktroub) call icopy(nn-ktroub,imap(ktroub+1),1,
     &                                imap(ktroub),1)
         nn = nn - 1
         print *,' remco lowers nn in few ',nn,ktroub
         it = 1
         iq = 0
         itr = itroub(ntroub)
         do 323 j=1,nequi
            do 232 k=2,nex(iq+1)+1
               it = it + 1
               if (it.eq.itr) then
                  write(iwr,602) itr-1,
     &                         (iex(1,m),iex(k,m),m=iq+1,iq+iequi(j))
602               format(/,' i do not like         ',i4,':',
     *            (t29,i4,'   =>',i4),'   it will be forbidden')
                  do 212 l=1,iequi(j)
                     if (iex(k,iq+l).le.nscf) then
                        nforbi = nforbi + 1
                        iforbi(nforbi,1) = iex(1,iq+l)
                        iforbi(nforbi,2) = iex(k,iq+l)
                     end if
212               continue
                  goto 121
               end if
232         continue
            iq = iq + iequi(j)
323      continue
121      continue
c....
c....    find out the partner degree of freedom (for print only)
c....
         do 2121 i=1,ktroub-1
            call icopy(i-1,imap(1)  ,1,imapc(1),1)
            call icopy(ktroub-i-1,imap(i+1),1,imapc(i),1)
            imapc(ktroub-1) = itroub(ntroub)
            call submat(hs,ktroub-1,imapc,ihfile,ihbloc,nstruc,'s')
c....
            call schmidt(hs,q,ktroub-1,ktr)
            if (ktr.eq.0) then
               if (i.eq.1) then
                  write(iwr,6111)
6111              format(/,' this brillouin state is (part of) the',
     &                   ' wavefunction')
                  goto 2121
               end if
               it = 1
               iq = 0
               itr = imap(i)
               do 3233 j=1,nequi
                  do 2323 k=2,nex(iq+1)+1
                     it = it + 1
                     if (it.eq.itr) then
                        write(iwr,6022)
     &                     itr-1,(iex(1,m),iex(k,m),m=iq+1,iq+iequi(j))
6022                    format(/,' matching excitation(s)',i4,':',
     &                  (t29,i4,'   =>',i4))
                     end if
2323              continue
                  iq = iq + iequi(j)
3233           continue
            end if
2121     continue
         if (nn.lt.1) call vberr('davidson collapses in few')
         goto 12345
      end if
c...     get h-submatrix
7001   call submat(hs,nn,imap,ihfile,ihbloc,nstruc,'h')
c      write (iwr,'(A)'),' h-submatrix'
c      call tripri(hs,nn)
c...     orthogonalise h-matrix
      call mult11(hs,hso,x,nn,nn,q,qt)
c...     diagonalise hso using orthog vectors as "unit-matrix"
      call jacobt(hso,iky,nn,q,nn,e,1,2,thresj,x)
c     call prvc(q,nn,nn,e,'v','a')
c
      if (struc) then
c...     doing the small ci (i.e. optimising the structures)
c...     see the =evec= directive in scfin (default: ivec=1)
         if (ivec.lt.1.or.ivec.gt.maxsel) call vberr
     &   ('wrong use of the evec directive')
         if (ivec.gt.1) then
c...        interesting...
            write(iwr,77) ivec
77          format(/,' optimising the orbitals for eigenvector :',i2)
            call prvc(q,nn,nn,e,'v','a')
         end if
         noffset=(ivec-1)*nn
      else
c...     doing the big ci (i.e. optimising the orbitals)
         noffset=0
      end if
c
      if (iprind.gt.20) then
         write(iwr,88) nn,e(1)+potnuc
88       format(/' trial energy for',i4,' dimensional subspace:',e27.14)
         write(iwr,99) (imap(i),q(i),i=1,nn)
99       format(' trial ci vector :',/,(t3,9(i5,f9.5)))
      end if
c
c...  fill appropriate elements of c with lowest eigenvector
c...  (or a higher one if requested by the =evec= directive)
c
      call vclr(c,1,nstruc)
      do 90 i=1,nn
         c(imap(i)) = q(i+noffset)
90    continue
c
c.... need to pick up c and some other housekeeping if other nodes bypass diag
      if (nn.lt.1) call vberr('zero trial ci-vector, degeneracies ?')
c
      if (ntroub.gt.0) then
c...     throw away trouble causing csf(s) effectively
         call elimin(ntroub,itroub,ihfile,ihbloc,nstruc)
      end if
c
      return
      end
      subroutine hsc(z,y,c,nstruc,ifile,ibloc,type)
c
      implicit real*8  (a-h,o-z) , integer   (i-n)
      character*1 type
c
c.....
c.....h/s * c is returned (h/s is real hermitian => triangular storage)
c.... uses labelled h and s matrices now
c....  efficency needs to be checked (FOKKE)
c.....
      common/matrixt/hof(170),sof(170),cricr(170),no,iblock,
     +               icof(170),irof(170)
      dimension z(nstruc),y(nstruc),c(nstruc)
      logical oz,oy
      external fget
c
      oz = type.eq.'h'.or.type.eq.'b'
      oy = type.eq.'s'.or.type.eq.'b'
c
      if (oz) call vclr(z,1,nstruc)
      if (oy) call vclr(y,1,nstruc)
c
      call search(ibloc,ifile)
10    call fget(hof,nw,ifile)
      if (nw.eq.0) return
       call unpack(cricr(1),32,icof,170)
       call unpack(cricr(86),32,irof,170)
      if (oz) then
         do 21 i=1,no
            z(icof(i)) = z(icof(i)) + hof(i)*c(irof(i))
            if (icof(i).ne.irof(i))
     1      z(irof(i)) = z(irof(i)) + hof(i)*c(icof(i))
21       continue
      end if
      if (oy) then
         do 22 i=1,no
            y(icof(i)) = y(icof(i)) + sof(i)*c(irof(i))
            if (icof(i).ne.irof(i)) 
     1      y(irof(i)) = y(irof(i)) + sof(i)*c(icof(i))
22       continue
      end if
      go to 10
      end    
      subroutine prhs(core,nstruc,ifile,ibloc,type)
c
      implicit real*8  (a-h,o-z) , integer   (i-n)
      dimension core(*)
c
c.....
c..... print h or s-matrix to nstruc structures
c..... they are  labelled stored on file ifile at block iblocf
c..... we need t get it in as a whole (!!!)
c..... derived from hcs
c.....
      common/matrixt/hsof(170,2),cricr(170),no,iblock,
     +               icof(170),irof(170)
      character*1 type
c
      
      integer ird, iwr, ipu, main, ibl2m, idaf, ibl3d, num8
      integer lenbl, ntapes, iblkci, numlib, iblkpl
      integer numdis, ibldis, irdcpy
      common /iofile/ ird,iwr,ipu,main,ibl2m,idaf,ibl3d,num8,
     +                lenbl(4),ntapes(5),iblkci,numlib,iblkpl,
     +                numdis,ibldis,irdcpy
c
      logical oroot
      external fget
      itri(i,j) = max(i,j)*(max(i,j)-1)/2+min(i,j)
c
      nn = nstruc*(nstruc+1)/2
c
      if (type.eq.'h'.or.type.eq.'1') then
         ityp = 1
      else if (type.eq.'s'.or.type.eq.'2') then
         ityp = 2
      else
         call vberr(' strange type in prhs')
      end if
      call vclr(core,1,nn)
c
      call search(ibloc,ifile)
10    call fget(hsof,nw,ifile)
      if (nw.eq.0)  go to 30
       call unpack(cricr(1),32,icof,170)
       call unpack(cricr(86),32,irof,170)
      do 20 i=1,no
         if (irof(i).gt.nstruc.or.icof(i).gt.nstruc) go to 20
         core(itri(irof(i),icof(i))) = hsof(i,ityp)
20    continue
      go to 10
c
c.... print 
c
30    it = 0
      do 40 i=1,nstruc
         write(iwr,600)  i,(core(it+j),j=1,i)
600      format(1x,i3,(t5,5f22.16))
40    it = it + i
c
      end
      subroutine geths(core,nstruc,ifile,ibloc,type)
c
      implicit real*8  (a-h,o-z) , integer   (i-n)
      dimension core(*)
c
c.....
c..... get h or s-matrix to nstruc structures
c..... they are  labelled stored on file ifile at block iblocf
c..... we need to get it in as a whole (!!!)
c..... derived from hcs
c.....
      common/matrixt/hsof(170,2),cricr(170),no,iblock,
     +               icof(170),irof(170)
      character*1 type
      external fget
c
      if (type.eq.'h'.or.type.eq.'1') then
         ityp = 1
      else if (type.eq.'s'.or.type.eq.'2') then
         ityp = 2
      else
         call vberr(' strange type in geths')
      end if
      nn =  nstruc*(nstruc+1)/2
      call vclr(core,1,nn)
c
      call search(ibloc,ifile)
10    call fget(hsof,nw,ifile)
      if (nw.eq.0)  go to 30
       call unpack(cricr(1),32,icof,170)
       call unpack(cricr(86),32,irof,170)
      do 20 i=1,no
         if (irof(i).gt.nstruc.or.icof(i).gt.nstruc) go to 20
         core(irof(i)*(irof(i)-1)/2+icof(i)) = hsof(i,ityp)
20    continue
      go to 10
c
30    continue
c
      end
      subroutine submat(a,n,imap,ifile,iblock,nstruc,type)
c
      implicit real*8  (a-h,o-z) , integer   (i-n)
      character*1 type
c
c.....
c.....a is constructed from the triangular matrix, stored on ifile,
c.....according to imap.
c.....
c..... cf writh and bamilt ; they have to supply nwidth in subinf
c
      dimension a(*),imap(*)
      common/matrixt/hsof(170,2),cricr(170),no,iblockk,
     +               icof(170),irof(170)
      common/subinf/ nwidth
c
      
      integer ird, iwr, ipu, main, ibl2m, idaf, ibl3d, num8
      integer lenbl, ntapes, iblkci, numlib, iblkpl
      integer numdis, ibldis, irdcpy
      common /iofile/ ird,iwr,ipu,main,ibl2m,idaf,ibl3d,num8,
     +                lenbl(4),ntapes(5),iblkci,numlib,iblkpl,
     +                numdis,ibldis,irdcpy
c
      external fget
      itri(i,j) = max(i,j)*(max(i,j)-1)/2 + min(i,j)
c
      if (type.eq.'h'.or.type.eq.'1') then
         ityp = 1
      else if (type.eq.'s'.or.type.eq.'2') then
         ityp = 2
      else
         call vberr(' strange type in submat')
      end if
      call vclr(a,1,n*(n+1)/2)
      maxc = min(imap(n)+nwidth,nstruc)
c
      call search(iblock,ifile)
10    call fget(hsof,nw,ifile)
      if (nw.eq.0)  go to 30
       call unpack(cricr(1),32,icof,170)
       call unpack(cricr(86),32,irof,170)
      do 20 i=1,no
c        if (irof(i).gt.maxc) go to 30
         do 15 j=1,n
            if (imap(j).eq.irof(i)) then
               do 14 k=1,n
                  if (imap(k).eq.icof(i)) then
                     a(itri(j,k)) = hsof(i,ityp)
                     go to 20
                  end if
14             continue
            end if
15       continue
20    continue
      go to 10
c
30    continue
      return
      end
      subroutine vbci(v,ci,maxvec)
c
      implicit real*8  (a-h,o-z), integer   (i-n)
c
c
c...   driver for vbci (normal vb) calculations
c...   orbitals are in v, q is scratch and returns vb-ci-vector
c
      dimension v(*),ci(*)
c
c     8note* tractlt noo included in vbdens.m mkd1vb
      integer iprint,nbasis,ncol,nsa,ncurt,lenact,lenbas,nactiv
      integer num3,iblkq,ionsec,nbsort,isecdu,n2int,setintpos
      logical incsort,oprs
      integer nnsort,isort,ksort,maxsort
      real*8 scri
c
      common /tractlt/ scri,iprint,nbasis,ncol,nsa,ncurt,lenact,lenbas,
     &                 num3,iblkq,ionsec,nbsort,isecdu,n2int,nactiv,
     &                 incsort,nnsort,isort,ksort,maxsort,oprs
c
c...  incsort : indicates if incore sorting is on
c...  nnsort : # sortbuffers
c...  isort : sortbuffer done
c...  ksort : adres of sortbuffers in qq (igmem_alloc)
c...  maxsort : maximum sort buffer used
      common/bypass/ index4,ihmat,idavid
      logical ohmat
      integer nstruc,iprinv,ihfile,ihbloc,ipbloc,idumhs
      common/hsinfo/ohmat(3),nstruc,iprinv,ihfile,ihbloc,ipbloc,idumhs
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
c...            symcri = critor*100 remcri = 1.0d-10 crihyb = 1.0d-14
c...            cribri = 1.0d-8 crinor = critor*100  cribub = 1.0d-10
c...            criign = 1.0d-44 cripop 1.0d-10
      real*8 critor,crilow,cridep,symcri,remcri,crihyb,cribri,
     1     crinor,cribub,criign,cripop
      common/vbcri/ critor,crilow,cridep,symcri,remcri,crihyb,cribri,
     1              crinor,cribub,criign,cripop
c
c... for dynamic memory allocation debugging purpose, contains IGMEM_ vars.
      integer mxheap
      parameter(mxheap=100)
      integer mxcount
      parameter(mxcount=100)

      integer iqoff, iq_heap, itag_heap
      integer igmem_count
      integer igmem_print
      integer igmem_priority
      integer igamem_priority
      integer igmem_size
      integer igamem_size
      integer igmem_totsize
      integer igamem_totsize
      integer igmem_maxsize
      integer igamem_maxsize
      logical ogmem_debug, ogmem_alloced_all, ogmem_nanify
      integer numheap
      integer numcount
      common/gmemdata/iqoff,iq_heap(mxheap),
     &     itag_heap(mxheap), igmem_count, numheap, 
     &     igmem_priority(mxheap),   igamem_priority(mxheap),
     &     igmem_size(mxheap),       igamem_size(mxheap),
     &     igmem_totsize(0:mxcount), igamem_totsize,
     &     igmem_maxsize(0:mxcount), igamem_maxsize, 
     &     igmem_print, numcount,
     &     ogmem_debug, ogmem_alloced_all, ogmem_nanify

      character*16 zgmem_varid,zgamem_arrnam
      common/gmemtext/zgmem_varid(mxheap),zgamem_arrnam(mxheap)
 
      integer mxresreg, mxressect, iadr, nres, nresreg,isize
      parameter(mxresreg=1)
      parameter(mxressect=40)
      common/gmemresv/iadr(0:mxressect,mxresreg),
     &     nres(mxresreg),isize(mxresreg),nresreg
      integer IGMEM_QUIET
      integer IGMEM_NORMAL
      integer IGMEM_DEBUG
      parameter (IGMEM_QUIET =-12)
      parameter (IGMEM_NORMAL=-11)
      parameter (IGMEM_DEBUG =-10)
c
c     ci :
c     ortvirt causes the virtuals to be mutually orthogonal in a ci
c     calculation. if (provirt) the occupied orbitals are projected out
c     ivirt is the first virtual
c     scf :
c     canonicalise (1): canonicalise doubles, (2): make virtuals diagonalise op
c     the operator is  0: none, 1: h, 2: fock
c     localise : pipec localise
c     idempotent : make projection operator idempotent
c     aos : suggests to use ao's to mix in
c
      common /vbvirt/ivirt,ortvirt,provirt,
     1               canonicalise,idempotent,aos,nit_loc
      logical ortvirt,provirt,idempotent,aos
      integer canonicalise(2),ivirt,nit_loc
c
c...  frozen core definition for TURTLE
c
      real*8 dxyz
      integer ncore,mapcie,mapiee,ncore_i,nsplice
      parameter (nsplice=512)
      common/splice/dxyz(3),ncore,ncore_i,
     1              mapcie(nsplice),mapiee(nsplice)
c
c...   turtle parameters
c
      integer melec,ncorq,maxex,maxact,mbonds,maxato,
     &        maxopa,ncases,maxcon,mmsocc,mxorbvb,maxspvb
c
      parameter (melec=300)
      parameter (ncorq=2000000)
      parameter (maxex=20000,maxact=150,mbonds=150,
     &           maxato=750,maxopa=1500)
      parameter (ncases=35)
c      ... maxcon in vbcrestr ...
c92      parameter (maxcon = 100)
      parameter (maxcon = 3000)
c      ... msocc max. # singly occuppieds ; was always enough (2009) ...
      parameter (mmsocc = 30)
c..   maxspvb  maximum # spinfunctions per configuration 
c..   (42 for 10 singlys and the ultimate answer)
      parameter (maxspvb=42) 
c
c     mxorbvb may be as large as allowed by n_8_16 (32767)
c
      parameter (mxorbvb=1000)
      real*8  shiscf,criscf,evb,ebi,swave,percmax,perccur,brmax,
     &      dampvb,fixb0
      integer optcri,iprins,nitscf,maxscf,isecv,nitdiis,maxdiis,ntscf,
     &        ipcmt,maxscf_s
      logical perall,scfconv,p0core,unitary
      common /scftvb/dampvb,fixb0,shiscf,criscf,evb,ebi,swave,percmax,
     &               perccur,brmax,
     &               optcri,iprins,nitscf,maxscf,isecv,nitdiis,
     &               maxdiis,ntscf,ipcmt,maxscf_s,perall,scfconv,
     &               p0core,unitary
      common/vbpert/crivarpert(4),shiftp,fpdoc,fpuoc,
     &              nvarp,ifirsp,idelp,npert,
     &              ntype,itype(mxorbvb),
     &              minfock,maxfock,minpert,maxpert,minvar,maxvarp,
     &              minfockd,maxfockd,
     &              fckp,ovbper,fockdiagonal
      INTEGER nvarp,ifirsp,idelp,npert,ntype,itype
      INTEGER minfock,maxfock,minpert,maxpert,minvar,maxvarp
      INTEGER minfockd,maxfockd
      logical ovbper,fckp,fockdiagonal
      real*8 crivarpert,shiftp,fpdoc,fpuoc
c description:
c   ntype -> 1 - specify optimisation per brillioun structure
c            2 - specify optimisation per orbital
c            3 - specify optimisation per kind (type excitation)
c   itype ->
c
c          ntype = 1
c            itype(1..maxex) - optimisation type for each (brillouin) structure
c
c          ntype = 2
c            itype(1..mxorbvb) - optimisation type for each orbital
c            
c          ntype = 3 
c            itype(1) - optimisation type for doubly to variably occupied
c            itype(2) - optimisation type for doubly to unoccupied
c            itype(3) - optimisation type for variably to variably occupied
c            itype(4) - optimisation type for variably to unoccupied
c            
c          ntype = 4
c            itype(1..mxorbvb) - optimisation type for each orbital (to be
c                               specified by clsbril
c
c   contents of itype array:
c            itype(x)=1 for Super CI optimalisation
c            itype(x)=2 for perturbation theory optimalisation
c            itype(x)=4 perturbation theory  determines choice between Super and Pert
c
c   crivarpert ~: critaria for itype=4 case per class if applicable

c
c..   final mainfile/integral info
c
c
c     core :  frozen core energy for the whole job
c     core_i :  frozen core energy  in makeps0 
c
      real*8 core,pot_nuc,core_i
      integer nfil,nofil,iblvb,lblvb,iblhs
      common /ffile/ core,core_i,pot_nuc,nfil,
     1               nofil(20),iblvb(20),lblvb(20),iblhs
      common/energi/enatom,nstate
      common /vbtimess/ tvirtu,t4indx0,t4indx,tmatre0,tmatre,tdavid,
     1                  tlagr,tgtr,tfock, 
     2                  n4indx0,n4indx,nmatre0,nmatre,nfock
      real*8  tvirtu,t4indx0,t4indx,tmatre0,tmatre,tdavid,tlagr,tgtr,
     3      tfock
      integer n4indx0,n4indx,nmatre0,nmatre,nfock

c
      
      integer ird, iwr, ipu, main, ibl2m, idaf, ibl3d, num8
      integer lenbl, ntapes, iblkci, numlib, iblkpl
      integer numdis, ibldis, irdcpy
      common /iofile/ ird,iwr,ipu,main,ibl2m,idaf,ibl3d,num8,
     +                lenbl(4),ntapes(5),iblkci,numlib,iblkpl,
     +                numdis,ibldis,irdcpy
c
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
      integer maxlablen
      parameter (maxlablen=20)
      character*20 lab
c
c  change the label definitions (init_time_periodsm in util3 )
c  also
c
      integer TP_ENTIRE
      parameter(TP_ENTIRE=1)
c
c  scf 
c
      integer TP_SCF
      parameter (TP_SCF=TP_ENTIRE+1)
      integer TP_ORFOG
      parameter (TP_ORFOG=TP_SCF+1)
      integer TP_DHSTAR
      parameter (TP_DHSTAR=TP_ORFOG+1)
      integer TP_RDMAT
      parameter (TP_RDMAT=TP_DHSTAR+1)
      integer TP_DIIS
      parameter (TP_DIIS=TP_RDMAT+1)
      integer TP_DIAG
      parameter (TP_DIAG=TP_DIIS+1)
c
c  mp2 code
c
      integer TP_APRDMP2
      parameter (TP_APRDMP2=TP_DIAG+1)

      integer TP_DHSTAR_GOP
      parameter (TP_DHSTAR_GOP=TP_APRDMP2+1)

c spare
      integer TP_APRM1234
      parameter (TP_APRM1234=TP_DHSTAR_GOP+1)
c spare
      integer TP_APRQ34
      parameter (TP_APRQ34=TP_APRM1234+1)

      integer TP_APRQ1
      parameter (TP_APRQ1=TP_APRQ34+1)
      integer TP_APRQ2
      parameter (TP_APRQ2=TP_APRQ1+1)
      integer TP_APRQ2D
      parameter (TP_APRQ2D=TP_APRQ2+1)
      integer TP_APRQ34D
      parameter (TP_APRQ34D=TP_APRQ2D+1)
      integer TP_APRMP2E
      parameter (TP_APRMP2E=TP_APRQ34D+1)

      integer TP_MP1PDM
      parameter (TP_MP1PDM=TP_APRMP2E+1)
      integer TP_MP1PDM_1
      parameter (TP_MP1PDM_1=TP_MP1PDM+1)
      integer TP_MP1PDM_2
      parameter (TP_MP1PDM_2=TP_MP1PDM_1+1)

      integer TP_APR1PDM
      parameter (TP_APR1PDM=TP_MP1PDM_2+1)

      integer TP_MP2HESS
      parameter (TP_MP2HESS=TP_APR1PDM+1)
      integer TP_MP2CHF
      parameter (TP_MP2CHF=TP_MP2HESS+1)
      integer TP_MP2MAKEW
      parameter (TP_MP2MAKEW=TP_MP2CHF+1)
      integer TP_MP2DS
      parameter (TP_MP2DS=TP_MP2MAKEW+1)
      integer TP_MP2BACK_P
      parameter (TP_MP2BACK_P=TP_MP2DS+1)
      integer TP_MP2BACK_F
      parameter (TP_MP2BACK_F=TP_MP2BACK_P+1)
      integer TP_MP2BACKTRAN_2
      parameter (TP_MP2BACKTRAN_2=TP_MP2BACK_F+1)
      integer TP_MP2MCDAB
      parameter (TP_MP2MCDAB=TP_MP2BACKTRAN_2+1)
c
c - parallel functions
c
      integer TP_DGOP
      parameter(TP_DGOP=TP_MP2MCDAB+1)

      integer TP_BCAST
      parameter(TP_BCAST=TP_DGOP+1)

      integer TP_NXTVAL
      parameter(TP_NXTVAL=TP_BCAST+1)

      integer TP_GENRAL
      parameter (TP_GENRAL=TP_NXTVAL+1)

      integer TP_GENRAL_1PDM
      parameter (TP_GENRAL_1PDM=TP_GENRAL+1)

      integer TP_GA_PUT_Q2D
      parameter (TP_GA_PUT_Q2D=TP_GENRAL_1PDM+1)

      integer TP_GA_ACC_Q2D
      parameter (TP_GA_ACC_Q2D=TP_GA_PUT_Q2D+1)

      integer TP_MKT2AO
      parameter (TP_MKT2AO   =TP_GA_ACC_Q2D+1)

      integer TP_APRL2
      parameter (TP_APRL2    =TP_MKT2AO+1)

      integer TP_GA_GET_L2
      parameter (TP_GA_GET_L2=TP_APRL2+1)

      integer TP_APRL34
      parameter (TP_APRL34   =TP_GA_GET_L2+1)

      integer TP_DGENRL
      parameter (TP_DGENRL   =TP_APRL34+1)

      integer TP_MP2
      parameter(TP_MP2       =TP_DGENRL+1)

      integer TP_JKDER
      parameter(TP_JKDER     =TP_MP2+1)

      integer TP_APRL1234
      parameter (TP_APRL1234 =TP_JKDER+1)

      integer TP_APRL1
      parameter (TP_APRL1    =TP_APRL1234+1)

      integer TP_APRDMP2_I
      parameter(TP_APRDMP2_I =TP_APRL1+1)

      integer TP_JKDER_GET
      parameter(TP_JKDER_GET =TP_APRDMP2_I+1)

      integer TP_MP1PDM_3
      parameter(TP_MP1PDM_3  =TP_JKDER_GET+1)

      integer TP_MP1PDM_4
      parameter(TP_MP1PDM_4  =TP_MP1PDM_3+1)

cgdf:  added TP_MKT2MO 14.3.95
      integer TP_MKT2MO
      parameter (TP_MKT2MO   =TP_MP1PDM_4+1)

csecd - time periods for second derivatives
c
      integer TP_2D_AOINTS
      parameter(TP_2D_AOINTS     =TP_MKT2MO+1)

      integer TP_2D_SCF
      parameter(TP_2D_SCF        =TP_2D_AOINTS+1)

      integer TP_2D_HFGRDN
      parameter(TP_2D_HFGRDN     =TP_2D_SCF+1)

      integer TP_2D_INDX2T
      parameter(TP_2D_INDX2T     =TP_2D_HFGRDN+1)

      integer TP_2D_MOINTS
      parameter(TP_2D_MOINTS     =TP_2D_INDX2T+1)

      integer TP_2D_TRNFKD
      parameter(TP_2D_TRNFKD     =TP_2D_MOINTS+1)

      integer TP_2D_CHFNDR
      parameter(TP_2D_CHFNDR     =TP_2D_TRNFKD+1)

      integer TP_2D_QMDER
      parameter(TP_2D_QMDER      =TP_2D_CHFNDR+1)

      integer TP_2D_DMDER
      parameter(TP_2D_DMDER      =TP_2D_QMDER+1)

      integer TP_2D_2D
      parameter(TP_2D_2D         =TP_2D_DMDER+1)

      integer TP_2D_CHF
      parameter(TP_2D_CHF        =TP_2D_2D+1)

      integer TP_2D_NUC
      parameter(TP_2D_NUC        =TP_2D_CHF+1)

      integer TP_2D_OVL
      parameter(TP_2D_OVL        =TP_2D_NUC+1)

      integer TP_2D_KE
      parameter(TP_2D_KE         =TP_2D_OVL+1)

      integer TP_2D_PE
      parameter(TP_2D_PE         =TP_2D_KE+1)

      integer TP_2D_2E
      parameter(TP_2D_2E         =TP_2D_PE+1)

      integer TP_2D_TOTAL
      parameter (TP_2D_TOTAL     =TP_2D_2E+1)

      integer TP_2D_CHFDRV
      parameter (TP_2D_CHFDRV    =TP_2D_TOTAL+1)

      integer TP_2D_PDENS
      parameter (TP_2D_PDENS     =TP_2D_CHFDRV+1)

      integer TP_2D_PFOCK
      parameter (TP_2D_PFOCK     =TP_2D_PDENS+1)

      integer TP_2D_CHFHESS
      parameter (TP_2D_CHFHESS   =TP_2D_PFOCK+1)

      integer TP_2D_CHFRHS
      parameter (TP_2D_CHFRHS    =TP_2D_CHFHESS+1)

      integer TP_2D_SYMMRHS
      parameter (TP_2D_SYMMRHS   =TP_2D_CHFRHS+1)

      integer TP_2D_SYMMU
      parameter (TP_2D_SYMMU     =TP_2D_SYMMRHS+1)

      integer TP_2D_PFOCK_OOOO
      parameter (TP_2D_PFOCK_OOOO=TP_2D_SYMMU+1)

      integer TP_2D_PFOCK_VOOO
      parameter (TP_2D_PFOCK_VOOO=TP_2D_PFOCK_OOOO+1)

      integer TP_2D_PFOCK_VVOO
      parameter (TP_2D_PFOCK_VVOO=TP_2D_PFOCK_VOOO+1)

      integer TP_2D_PFOCK_VOVO
      parameter (TP_2D_PFOCK_VOVO=TP_2D_PFOCK_VVOO+1)

      integer TP_2D_PFOCK_SUM
      parameter (TP_2D_PFOCK_SUM =TP_2D_PFOCK_VOVO+1)

      integer TP_2D_AOGEN
      parameter(TP_2D_AOGEN=TP_2D_PFOCK_SUM+1)

      integer TP_2D_AOUT
      parameter(TP_2D_AOUT =TP_2D_AOGEN+1)

      integer TP_TEST1
      parameter(TP_TEST1   =TP_2D_AOUT+1)

      integer TP_TEST2
      parameter(TP_TEST2   =TP_TEST1+1)

      integer TP_TEST3
      parameter(TP_TEST3   =TP_TEST2+1)

      integer TP_TEST4
      parameter(TP_TEST4   =TP_TEST3+1)

      integer TP_TEST5
      parameter(TP_TEST5   =TP_TEST4+1)

      integer TP_TEST6
      parameter(TP_TEST6   =TP_TEST5+1)

      integer TP_HFGRAD
      parameter(TP_HFGRAD  =TP_TEST6+1)

      integer TP_GAMULT2
      parameter(TP_GAMULT2 =TP_HFGRAD+1)

      integer TP_GAORTHOG
      parameter(TP_GAORTHOG=TP_GAMULT2+1)

      integer TP_PDIAG
      parameter (TP_PDIAG  =TP_GAORTHOG+1)

      integer TP_MULT2
      parameter (TP_MULT2  =TP_PDIAG+1)

      integer TP_INTEG
      parameter (TP_INTEG  =TP_MULT2+1)
c
c =================  I/O timers ========================
c
c find
c	
      integer TP_IOFM1, TP_IOF0, TP_IOF1, TP_IOF2, TP_IOF3,
     &	TP_IOF4, TP_IOF5, TP_IOF6, TP_IOF7

      parameter(TP_IOFM1=TP_INTEG+1)
      parameter (TP_IOF0=TP_IOFM1+1)
      parameter (TP_IOF1=TP_IOF0+1)
      parameter (TP_IOF2=TP_IOF1+1)
      parameter (TP_IOF3=TP_IOF2+1)
      parameter (TP_IOF4=TP_IOF3+1)
      parameter (TP_IOF5=TP_IOF4+1)
      parameter (TP_IOF6=TP_IOF5+1)
      parameter (TP_IOF7=TP_IOF6+1)
c
c get
c
      integer TP_IOGM1, TP_IOG0, TP_IOG1, TP_IOG2, TP_IOG3, 
     &  TP_IOG4, TP_IOG5, TP_IOG6, TP_IOG7
      parameter (TP_IOGM1=TP_IOF7+1)
      parameter (TP_IOG0=TP_IOGM1+1)
      parameter (TP_IOG1=TP_IOG0+1)
      parameter (TP_IOG2=TP_IOG1+1)
      parameter (TP_IOG3=TP_IOG2+1)
      parameter (TP_IOG4=TP_IOG3+1)
      parameter (TP_IOG5=TP_IOG4+1)
      parameter (TP_IOG6=TP_IOG5+1)
      parameter (TP_IOG7=TP_IOG6+1)
c
c put
c
      integer TP_IOPM1,  TP_IOP0, TP_IOP1, TP_IOP2, TP_IOP3,
     & TP_IOP4, TP_IOP5, TP_IOP6, TP_IOP7
      parameter (TP_IOPM1=TP_IOG7+1)
      parameter (TP_IOP0=TP_IOPM1+1)
      parameter (TP_IOP1=TP_IOP0+1)
      parameter (TP_IOP2=TP_IOP1+1)
      parameter (TP_IOP3=TP_IOP2+1)
      parameter (TP_IOP4=TP_IOP3+1)
      parameter (TP_IOP5=TP_IOP4+1)
      parameter (TP_IOP6=TP_IOP5+1)
      parameter (TP_IOP7=TP_IOP6+1)
c
c open
c
      integer TP_IOOM1,TP_IOO0,TP_IOO1,TP_IOO2,TP_IOO3,
     & TP_IOO4,TP_IOO5, TP_IOO6, TP_IOO7

      parameter (TP_IOOM1=TP_IOP7+1)
      parameter (TP_IOO0=TP_IOOM1+1)
      parameter (TP_IOO1=TP_IOO0+1)
      parameter (TP_IOO2=TP_IOO1+1)
      parameter (TP_IOO3=TP_IOO2+1)
      parameter (TP_IOO4=TP_IOO3+1)
      parameter (TP_IOO5=TP_IOO4+1)
      parameter (TP_IOO6=TP_IOO5+1)
      parameter (TP_IOO7=TP_IOO6+1)
c
c delfil (only significant for GA-files dumped to disc
c
      integer TP_IO_GAFILE_READ, TP_IO_GAFILE_DUMP
      parameter (TP_IO_GAFILE_READ      =TP_IOO7+1)
      parameter (TP_IO_GAFILE_DUMP      =TP_IO_GAFILE_READ+1)
c
c Peigs parallel diag
c
      integer TP_PEIGS
      parameter (TP_PEIGS               =TP_IO_GAFILE_DUMP+1)
c
c Scalapack parallel diag
c
      integer TP_PDSYEV
      parameter (TP_PDSYEV               =TP_PEIGS+1)
      integer TP_PDSYEVX
      parameter (TP_PDSYEVX              =TP_PDSYEV+1)
      integer TP_PDSYEVD
      parameter (TP_PDSYEVD              =TP_PDSYEVX+1)
      integer TP_PDSYEVR
      parameter (TP_PDSYEVR              =TP_PDSYEVD+1)
c
c timers for CCP1 DFT module
c
      integer    TP_DFT_JFIT
      parameter (TP_DFT_JFIT            =TP_PDSYEVR+1)
      integer    TP_DFT_JFIT_VFORM
      parameter (TP_DFT_JFIT_VFORM      =TP_DFT_JFIT+1)
      integer    TP_DFT_JFIT_TR
      parameter (TP_DFT_JFIT_TR         =TP_DFT_JFIT_VFORM+1)
      integer    TP_DFT_JFIT_NR
      parameter (TP_DFT_JFIT_NR         =TP_DFT_JFIT_TR+1)
      integer    TP_DFT_JFIT_COEF
      parameter (TP_DFT_JFIT_COEF       =TP_DFT_JFIT_NR+1)
      integer    TP_DFT_JFIT_KSMAT
      parameter (TP_DFT_JFIT_KSMAT      =TP_DFT_JFIT_COEF+1)
      integer    TP_DFT_JFIT_ENERGY
      parameter (TP_DFT_JFIT_ENERGY     =TP_DFT_JFIT_KSMAT+1)
      integer    TP_DFT_EXQUAD
      parameter (TP_DFT_EXQUAD          =TP_DFT_JFIT_ENERGY+1)
      integer    TP_DFT_EXQUAD_INTRO
      parameter (TP_DFT_EXQUAD_INTRO    =TP_DFT_EXQUAD+1)
      integer    TP_DFT_EXQUAD_INTEG
      parameter (TP_DFT_EXQUAD_INTEG    =TP_DFT_EXQUAD_INTRO+1)
      integer    TP_DFT_EXQUAD_DGOP
      parameter (TP_DFT_EXQUAD_DGOP     =TP_DFT_EXQUAD_INTEG+1)
      integer    TP_DFT_EXQUADF
      parameter (TP_DFT_EXQUADF         =TP_DFT_EXQUAD_DGOP+1)
      integer    TP_DFT_EXQUADF_INTRO
      parameter (TP_DFT_EXQUADF_INTRO   =TP_DFT_EXQUADF+1)
      integer    TP_DFT_EXQUADF_INTEG
      parameter (TP_DFT_EXQUADF_INTEG   =TP_DFT_EXQUADF_INTRO+1)
      integer    TP_DFT_EXQUADF_DGOP
      parameter (TP_DFT_EXQUADF_DGOP    =TP_DFT_EXQUADF_INTEG+1)
      integer    TP_DFT_EXQUADLHS
      parameter (TP_DFT_EXQUADLHS       =TP_DFT_EXQUADF_DGOP+1)
      integer    TP_DFT_EXQUADLHS_INTRO
      parameter (TP_DFT_EXQUADLHS_INTRO =TP_DFT_EXQUADLHS+1)
      integer    TP_DFT_EXQUADLHS_INTEG
      parameter (TP_DFT_EXQUADLHS_INTEG =TP_DFT_EXQUADLHS_INTRO+1)
      integer    TP_DFT_EXQUADLHS_DGOP
      parameter (TP_DFT_EXQUADLHS_DGOP  =TP_DFT_EXQUADLHS_INTEG+1)
      integer    TP_DFT_EXQUADHES
      parameter (TP_DFT_EXQUADHES       =TP_DFT_EXQUADLHS_DGOP+1)
      integer    TP_DFT_EXQUADHES_INTRO
      parameter (TP_DFT_EXQUADHES_INTRO =TP_DFT_EXQUADHES+1)
      integer    TP_DFT_EXQUADHES_INTEG
      parameter (TP_DFT_EXQUADHES_INTEG =TP_DFT_EXQUADHES_INTRO+1)
      integer    TP_DFT_EXQUADHES_DGOP
      parameter (TP_DFT_EXQUADHES_DGOP  =TP_DFT_EXQUADHES_INTEG+1)
      integer    TP_DFT_EXQUADRHS
      parameter (TP_DFT_EXQUADRHS       =TP_DFT_EXQUADHES_DGOP+1)
      integer    TP_DFT_EXQUADRHS_INTRO
      parameter (TP_DFT_EXQUADRHS_INTRO =TP_DFT_EXQUADRHS+1)
      integer    TP_DFT_EXQUADRHS_INTEG
      parameter (TP_DFT_EXQUADRHS_INTEG =TP_DFT_EXQUADRHS_INTRO+1)
      integer    TP_DFT_EXQUADRHS_DGOP
      parameter (TP_DFT_EXQUADRHS_DGOP  =TP_DFT_EXQUADRHS_INTEG+1)
      integer    TP_DFT_EXQUADDKSX
      parameter (TP_DFT_EXQUADDKSX      =TP_DFT_EXQUADRHS_DGOP+1)
      integer    TP_DFT_EXQUADDKSX_INTRO
      parameter (TP_DFT_EXQUADDKSX_INTRO=TP_DFT_EXQUADDKSX+1)
      integer    TP_DFT_EXQUADDKSX_INTEG
      parameter (TP_DFT_EXQUADDKSX_INTEG=TP_DFT_EXQUADDKSX_INTRO+1)
      integer    TP_DFT_EXQUADDKSX_DGOP
      parameter (TP_DFT_EXQUADDKSX_DGOP =TP_DFT_EXQUADDKSX_INTEG+1)
      integer    TP_DFT_EXQUADDKS
      parameter (TP_DFT_EXQUADDKS       =TP_DFT_EXQUADDKSX_DGOP+1)
      integer    TP_DFT_EXQUADDKS_INTRO
      parameter (TP_DFT_EXQUADDKS_INTRO =TP_DFT_EXQUADDKS+1)
      integer    TP_DFT_EXQUADDKS_INTEG
      parameter (TP_DFT_EXQUADDKS_INTEG =TP_DFT_EXQUADDKS_INTRO+1)
      integer    TP_DFT_EXQUADDKS_DGOP
      parameter (TP_DFT_EXQUADDKS_DGOP  =TP_DFT_EXQUADDKS_INTEG+1)



      integer    TP_DFT_JMULT
      parameter (TP_DFT_JMULT           =TP_DFT_EXQUADDKS_DGOP+1)
      integer    TP_DFT_JMULT_INTRO
      parameter (TP_DFT_JMULT_INTRO     =TP_DFT_JMULT+1)
      integer    TP_DFT_JMULT_SB
      parameter (TP_DFT_JMULT_SB        =TP_DFT_JMULT_INTRO+1)
      integer    TP_DFT_JMULT_FOCK
      parameter (TP_DFT_JMULT_FOCK      =TP_DFT_JMULT_SB+1)
      integer    TP_DFT_JMULT_CJAT0
      parameter (TP_DFT_JMULT_CJAT0     =TP_DFT_JMULT_FOCK+1)
c
c coulomb fitted gradients
c
      integer    TP_DFT_JFITG
      parameter (TP_DFT_JFITG           =TP_DFT_JMULT_CJAT0+1)

      integer    TP_DFT_JFITG_VFORM
      parameter (TP_DFT_JFITG_VFORM     =TP_DFT_JFITG+1)
      integer    TP_DFT_JFITG_TR
      parameter (TP_DFT_JFITG_TR        =TP_DFT_JFITG_VFORM+1)
      integer    TP_DFT_JFITG_NR
      parameter (TP_DFT_JFITG_NR        =TP_DFT_JFITG_TR+1)
      integer    TP_DFT_JFITG_COEF
      parameter (TP_DFT_JFITG_COEF      =TP_DFT_JFITG_NR+1)
      integer    TP_DFT_JFITG_2C
      parameter (TP_DFT_JFITG_2C        =TP_DFT_JFITG_COEF+1)
      integer    TP_DFT_JFITG_3C
      parameter (TP_DFT_JFITG_3C        =TP_DFT_JFITG_2C+1)
      integer    TP_DFT_JFIT_TR_INIT
      parameter (TP_DFT_JFIT_TR_INIT    =TP_DFT_JFITG_3C+1)
      integer    TP_DFT_JFIT_INV
      parameter (TP_DFT_JFIT_INV        =TP_DFT_JFIT_TR_INIT+1)
c
c VB
c
      integer TP_VB
      parameter (TP_VB                  =TP_DFT_JFIT_INV+1)
      integer TP_VB_STRUC
      parameter (TP_VB_STRUC            =TP_VB+1)
      integer TP_VB_ME
      parameter (TP_VB_ME               =TP_VB_STRUC+1)
      integer TP_VB_DIAG
      parameter (TP_VB_DIAG             =TP_VB_ME+1)
      integer TP_VB_TRAN
      parameter (TP_VB_TRAN             =TP_VB_DIAG+1)
      integer TP_VB_VIRT
      parameter (TP_VB_VIRT             =TP_VB_TRAN+1)
      integer TP_VB_LADM
      parameter (TP_VB_LADM             =TP_VB_VIRT+1)
      integer TP_VB_DTRAN
      parameter (TP_VB_DTRAN            =TP_VB_LADM+1)
c
c VB parallel extra
c
      integer TP_ISEND
      parameter (TP_ISEND               =TP_VB_DTRAN+1)
      integer TP_IRECV
      parameter (TP_IRECV               =TP_ISEND+1)
      integer TP_WAIT
      parameter (TP_WAIT                =TP_IRECV+1)
c
c One electron derivatives
c
      integer TP_STVECP
      parameter (TP_STVECP              =TP_WAIT+1)

      integer TP_TVDER
      parameter (TP_TVDER               =TP_STVECP+1)

      integer TP_SDER
      parameter (TP_SDER                =TP_TVDER+1)

      integer TP_SGRAD
      parameter (TP_SGRAD               =TP_SDER+1)

      integer TP_HELFEY
      parameter (TP_HELFEY              =TP_SGRAD+1)


      !F90 time periods start here
      integer TP_F90_START
      parameter( TP_F90_START           = TP_HELFEY+1 )
      integer TP_F90_SCF
      parameter( TP_F90_SCF             = TP_F90_START+1 )
      integer TP_F90_BUILD
      parameter( TP_F90_BUILD           = TP_F90_SCF+1 )
      integer TP_F90_DIIS
      parameter( TP_F90_DIIS            = TP_F90_BUILD+1 )
      integer TP_F90_SIMIL
      parameter( TP_F90_SIMIL           = TP_F90_DIIS+1 )
      integer TP_F90_DIAG
      parameter( TP_F90_DIAG            = TP_F90_SIMIL+1 )
      integer TP_F90_BACK
      parameter( TP_F90_BACK            = TP_F90_DIAG+1 )
      integer TP_F90_ASSIGN
      parameter( TP_F90_ASSIGN          = TP_F90_BACK+1 )
      integer TP_F90_ORTHOG
      parameter( TP_F90_ORTHOG          = TP_F90_ASSIGN+1 )
      integer TP_F90_MAKE_DENS
      parameter( TP_F90_MAKE_DENS       = TP_F90_ORTHOG+1 )
      integer TP_F90_END
      parameter( TP_F90_END             = TP_F90_MAKE_DENS+1 )
      integer TP_F90_LEV_SHIFT
      parameter( TP_F90_LEV_SHIFT       = TP_F90_END+1 )
      integer TP_F90_TESTER_EVAL
      parameter( TP_F90_TESTER_EVAL     = TP_F90_LEV_SHIFT+1 )
      integer TP_F90_DELTA_EVAL
      parameter( TP_F90_DELTA_EVAL      = TP_F90_TESTER_EVAL+1 )
      integer TP_F90_TDOWN
      parameter( TP_F90_TDOWN           = TP_F90_DELTA_EVAL+1 )
      integer TP_F90_RDMAT
      parameter( TP_F90_RDMAT           = TP_F90_TDOWN+1 )
      integer TP_F90_INTS
      parameter( TP_F90_INTS            = TP_F90_RDMAT+1 )
      integer TP_NEWSCF
      parameter( TP_NEWSCF              = TP_F90_INTS+1 )

      integer TP_DENSCF
      parameter( TP_DENSCF              = TP_NEWSCF+1 )
      integer TP_DENSCF_BUILD
      parameter( TP_DENSCF_BUILD        = TP_DENSCF+1 )
      integer TP_DENSCF_RDMAT
      parameter( TP_DENSCF_RDMAT        = TP_DENSCF_BUILD+1 )
      integer TP_DENSCF_INTS
      parameter( TP_DENSCF_INTS         = TP_DENSCF_RDMAT+1 )
      integer TP_DENSCF_DIAG_S
      parameter( TP_DENSCF_DIAG_S       = TP_DENSCF_INTS+1 )
      integer TP_DENSCF_SIMIL
      parameter( TP_DENSCF_SIMIL        = TP_DENSCF_DIAG_S+1 )
      integer TP_DENSCF_DIAG
      parameter( TP_DENSCF_DIAG         = TP_DENSCF_SIMIL+1 )
      integer TP_DENSCF_BACK
      parameter( TP_DENSCF_BACK         = TP_DENSCF_DIAG+1 )
      integer TP_DENSCF_MAKE_DENS
      parameter( TP_DENSCF_MAKE_DENS    = TP_DENSCF_BACK+1 )
      integer TP_DENSCF_TDOWN
      parameter( TP_DENSCF_TDOWN        = TP_DENSCF_MAKE_DENS+1 )

      integer TP_DRHFCL_GA
      parameter( TP_DRHFCL_GA           = TP_DENSCF_TDOWN+1 )
c
c     RPA module
c
      integer TP_RESPONSE
      parameter( TP_RESPONSE            = TP_DRHFCL_GA + 1)
      integer TP_RPA
      parameter( TP_RPA                 = TP_RESPONSE + 1)
      integer TP_TDA
      parameter( TP_TDA                 = TP_RPA + 1)
      integer TP_RPANAL
      parameter( TP_RPANAL              = TP_TDA + 1)
      integer TP_RPA_MO2AO
      parameter( TP_RPA_MO2AO           = TP_RPANAL + 1)
      integer TP_RPA_INT
      parameter( TP_RPA_INT             = TP_RPA_MO2AO + 1)
      integer TP_RPA_CNTRCT
      parameter( TP_RPA_CNTRCT          = TP_RPA_INT + 1)
      integer TP_RPA_AO2MO
      parameter( TP_RPA_AO2MO           = TP_RPA_CNTRCT + 1)
      integer TP_TDA_MO2AO
      parameter( TP_TDA_MO2AO           = TP_RPA_AO2MO + 1)
      integer TP_TDA_INT
      parameter( TP_TDA_INT             = TP_TDA_MO2AO + 1)
      integer TP_TDA_CNTRCT
      parameter( TP_TDA_CNTRCT          = TP_TDA_INT + 1)
      integer TP_TDA_AO2MO
      parameter( TP_TDA_AO2MO           = TP_TDA_CNTRCT + 1)
c
c     Define the common blocks
c
      integer maxtp
      parameter (maxtp=TP_TDA_AO2MO)
      integer maxtpdepth
      parameter (maxtpdepth = 10)
      integer itpdepth
      integer itpstack
      integer ntpc
      integer parent
      real*8  ttotw, ttotc, tsw, tsc
      real*8  ttotu, ttots, tsu, tss
      real*8  taggc
      common/timeperiods/ttotw(maxtp),ttotc(maxtp),
     &     ttotu(maxtp),ttots(maxtp),
     &     tsw(maxtp),tsc(maxtp),
     &     tsu(maxtp),tss(maxtp),
     &     taggc(maxtp),
     &     ntpc(maxtp),parent(maxtp),
     &     itpstack(0:maxtpdepth),itpdepth
      common/timeperiodsc/lab(maxtp)
      real*8 qq
      integer ivoff, memhandle
      common/vcore/qq(2)
      common/vcoreoff/ivoff, memhandle
      logical struc
      common/hcpu_vb/iprhs,irowhs
      ind(i,j) = max(i,j)*(max(i,j)-1)/2 + min(i,j)
      struc = .true. 
      ivec = 1
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
      call getqvb(v,nbasis,ncol,isecv,'nopr')
c
c...    reserve space on ed7 for vectors  (dumped in  transformvb)
c
      kvb7_transvb = kscra7vb('kvb7_transvb',nbasis*ncol,'r','w')
c
         call wr15vb(ndum1,ndum2,ndum3,ndum4,ndum5,ndum6,ndum7,
     &               ndum8,ndum9,imax,ndum11,ndum12,'read')

      call vbfirst(v,maxvec)
c
c...  4-index-transformation
c
      a1 = cpulft(1)
      call start_time_period(TP_VB_TRAN)

      ng = lenact*(lenact+1)/2  + 1
      ks = igmem_alloc_inf(lenact,'vbci.m','vbci','ks',IGMEM_DEBUG)
      kh = igmem_alloc_inf(lenact,'vbci.m','vbci','kh',IGMEM_DEBUG)
c
c...  (first integral = (00/00) = 0.0 )

      if (index4.eq.1) then
         if (provirt.or.ortvirt) then
c...        do some orthogonalisations
            nnot = 0
            if (.not.provirt) nnot = ivirt - 1
            ksao = igmem_alloc_inf(nbasis*(nbasis+1)/2,'vbci.m','vbci',
     &                             'ksao',IGMEM_DEBUG)
            kvcopy = igmem_alloc_inf(nbasis*ncol,'vbci.m','vbci',
     &                               'kvcopy',IGMEM_DEBUG)
            kmaxmem = igmem_max_memory()
            kmemscr = kmaxmem/1000
            kscr = igmem_alloc_inf(kmemscr,'vbci.m','vbci','kscr',
     &                             IGMEM_DEBUG)
            call get1e(qq(ksao),dummy,'s',qq(kscr))
            nvect = ncol - ncore - nnot
            nvec  = ncol - ncore - ivirt + 1
            call fmove(v((ncore+nnot)*nbasis+1),qq(kvcopy),nvect*nbasis)
            call normvc(qq(kvcopy),qq(ksao),qq(kscr),nbasis,nvect,
     &                  crilow)
            call gmem_free_inf(kscr,'vbci.m','vbci','kscr')
            if (nnot.eq.0) then
               call fmove(qq(kvcopy+(ivirt-1)*nbasis),
     &                    v((ncore+ivirt-1)*nbasis+1),nvec*nbasis)
            else
               call fmove(qq(kvcopy),
     &                    v((ncore+ivirt-1)*nbasis+1),nvec*nbasis)
            end if
            call gmem_free_inf(kvcopy,'vbci.m','vbci','kvcopy')
            call gmem_free_inf(ksao,'vbci.m','vbci','ksao')
         end if
         call transformvb(qq(ks),qq(kh),v)
      end if
      kg = igmem_alloc_inf(ng,'vbci.m','vbci','kg',IGMEM_DEBUG)
      call getin2(qq(kg))
      call clredx
c
c...  restore integrals if needed
c
      if (ihmat.eq.1.and.index4.ne.1) then
         call getin2(qq(kg))
         call rdedx(qq(kh),lenact,iblhs,num3)
         call reads(qq(ks),lenact,num3)
      end if
c
c...  print integrals if requested  (iprint)
c
      a1 = cpulft(1) - a1
      t4indx = a1
      n4indx = 1
      call end_time_period(TP_VB_TRAN)
      if (iprint.ge.50) then
601      format(/'  ---',a45,'---',/)
         write(iwr,601) ' 1-electron matrix between active orbitals '
         call tripri(qq(kh),nsa)
         write(iwr,601) ' overlap matrix between active orbitals '
         call tripri(qq(ks),nsa)
         if (iprint.ge.1000000) then
            write(iwr,601) ' 2-electron matrix between active orbitals '
            call tripri(qq(kg+1),lenact)
         end if
      end if

c
c...   build h/s-matrix
c
      a1 = cpulft(1)
      call start_time_period(TP_VB_ME)

      if (ihmat.eq.1) then
         kmaxmem = igmem_max_memory()
         kmemscr = kmaxmem/10

         call hsmat(qq(ks),qq(kh),qq(kg),nsa,qq)
      end if
      call gmem_free_inf(kg,'vbci.m','vbci','kg')
      call gmem_free_inf(kh,'vbci.m','vbci','kh')
      call gmem_free_inf(ks,'vbci.m','vbci','ks')

      a1 = cpulft(1) - a1
      tmatre = a1
      nmatre = 1
      call end_time_period(TP_VB_ME)
c
c...  f only partial matrix produced
c
      if (irowhs.gt.1) return
c...
c...  diagonalise h-matrix
c
      a1 = cpulft(1)
      call start_time_period(TP_VB_DIAG)
      if (idavid.eq.1) then
         kmaxmem = igmem_max_memory()
         kmemscr = kmaxmem/10
         kscr = igmem_alloc_inf(kmemscr,'vbci.m','vbci','kscr',
     &                          IGMEM_DEBUG)

         call davidt(qq(kscr),kmemscr,core,evb,'vbci')
         do i=1,nstruc
           ci(i) = qq(kscr+i-1)
         end do
         call gmem_free_inf(kscr,'vbci.m','vbci','kscr')
      end if
c
      write(iwr,76) evb
76    format(1x,/,' ==================================================',
     1          /,'    VB-energy from VBCI',f24.14, 
     2          /,' ==================================================')
      a1 = cpulft(1) - a1
      tdavid = a1
      call end_time_period(TP_VB_DIAG)
      if (nstruc.le.maxcon) then

         ki = igmem_alloc_inf(nstruc,'vbci.m','vbci','ki',IGMEM_DEBUG)
         ks = igmem_alloc_inf(nstruc*(nstruc+1)/2,'vbci.m','vbci',
     &                        'ks',IGMEM_DEBUG)
         call filnat(qq(ki),nstruc)
         call submat(qq(ks),nstruc,qq(ki),ihfile,ihbloc,nstruc,'s')

         kwgn = igmem_alloc_inf(nstruc,'vbci.m','vbci','kwgn',
     &       IGMEM_DEBUG)
         ksscr = igmem_alloc_inf(nstruc*nstruc,'vbci.m','vbci',
     &        'ksscr',IGMEM_DEBUG)
         ksinv = igmem_alloc_inf(nstruc*nstruc,'vbci.m','vbci',
     &        'ksinv',IGMEM_DEBUG)
         ksv1 = igmem_alloc_inf(nstruc,'vbci.m','vbci','ksv1',
     &       IGMEM_DEBUG)
         ksv2 = igmem_alloc_inf(nstruc,'vbci.m','vbci','ksv2',
     &       IGMEM_DEBUG)

         call submat(qq(ksscr),nstruc,qq(ki),ihfile,ihbloc,nstruc,'s')
         call square(qq(ksinv),qq(ksscr),nstruc,nstruc)

         fcrit = 10E-8
         call osinv_vb(qq(ksinv),nstruc,idetval,fcrit,qq(ksv1),qq(ksv2))

         it = 1
         psinorm = 0.0d0
         do ij=1,nstruc
           psinorm = psinorm + ci(ij)*ci(ij)/qq(ksinv+it-1)
           it = it + nstruc + 1
         end do

         it = 1
         do ij=1,nstruc
           qq(kwgn+ij-1) = ci(ij)*ci(ij)/qq(ksinv+it-1)
           qq(kwgn+ij-1) = qq(kwgn+ij-1)/psinorm
           it = it + nstruc + 1
         end do

         kww = igmem_alloc_inf(nstruc,'vbci.m','vbci','kww',IGMEM_DEBUG)
         do 700 i=1,nstruc
            qq(kww+i-1) = 0.0d0
            do 699 j=1,nstruc
               qq(kww+i-1) = qq(kww+i-1) + ci(i)*qq(ks+ind(i,j)-1)*ci(j)
699         continue
700      continue
         write(iwr,607) (ci(ij),ij=1,nstruc)
607      format(/,' vb-vector :',(t13,10f10.6))
         write(iwr,771) (qq(kww +ij-1),ij=1,nstruc)
         write(iwr,772) (qq(kwgn+ij-1),ij=1,nstruc)
771      format(/,' weights CC:',(t13,10f10.6))
772      format(/,' weights GN:',(t13,10f10.6))

         call wr15vb(nelec,ndum2,nconf,ndum4,ndum5,ndum6,ndum7,
     &            ndum8,ndum9,ndum10,ndum11,ndum12,'read')
         write(iwr,78) (qq(ks+ind(i,i)-1),i=1,nstruc)
78       format(/,' struc-norm:',(t13,10f10.6))
      else 
         kww = igmem_alloc_inf(nstruc,'vbci.m','vbci','kww',IGMEM_DEBUG)
      end if

      kbonds = igmem_alloc_inf(nelec*nstruc,'vbci.m','vbci','kbonds',
     &                         IGMEM_DEBUG)
      kmemscr = 2*nstruc + 4*nstruc*nstruc + 
     &          2*(nstruc*(nstruc+1)/2) + 1 +
     &          nstruc*nelec
      if ( nstruc .gt. 20 ) then
        kmemscr = nstruc + nstruc*nstruc +
     &            2*(nstruc*(nstruc+1)/2) + 1 +
     &            nstruc*nelec
      end if
      kscr = igmem_alloc_inf(kmemscr,'vbci.m','vbci','kscr',
     &                       IGMEM_DEBUG)
      call prweig(qq(kww),qq(kbonds),qq(kscr),nstruc,qq(kwgn))


      kcin = igmem_alloc_inf(nstruc,'vbci.m','vbci','kcin',
     &       IGMEM_DEBUG)

      if (nstruc.le.20) nstate = max(nstruc,nstate)
      if (nstate.ne.0) call prhsmat(nstate,qq(kscr),evb,ci,qq(kcin))
      if (nstruc.gt.20) then
        call prheigval(nstate,qq(kscr),qq(kwgn))
      end if

      call gmem_free_inf(kcin,'vbci.m','vbci','kcin')
      call gmem_free_inf(kscr,'vbci.m','vbci','kscr')
      call gmem_free_inf(kbonds,'vbci.m','vbci','kbonds')
      call gmem_free_inf(kww,'vbci.m','vbci','kww')

      call gmem_free_inf(ksv2,'vbci.m','vbci','ksv2')
      call gmem_free_inf(ksv1,'vbci.m','vbci','ksv1')
      call gmem_free_inf(ksinv,'vbci.m','vbci','ksinv')
      call gmem_free_inf(ksscr,'vbci.m','vbci','ksscr')
      call gmem_free_inf(kwgn,'vbci.m','vbci','kwgn')

      call gmem_free_inf(ks,'vbci.m','vbci','ks')
      call gmem_free_inf(ki,'vbci.m','vbci','ki')
      call putqvb(v,nbasis,ncore+nsa)
c
c...  determine natural orbitals
c
      call buildnat(ci,qq)

      return
      end
      subroutine buildnat(ci,qq)
c
c...   interface for natorbt ;- mig6ht include call to psi0det sometime
c
      implicit real*8 (a-h,o-z), integer (i-n)
c
      dimension ci(*),qq(*)
c
      integer mxheap
      parameter(mxheap=100)
      integer mxcount
      parameter(mxcount=100)

      integer iqoff, iq_heap, itag_heap
      integer igmem_count
      integer igmem_print
      integer igmem_priority
      integer igamem_priority
      integer igmem_size
      integer igamem_size
      integer igmem_totsize
      integer igamem_totsize
      integer igmem_maxsize
      integer igamem_maxsize
      logical ogmem_debug, ogmem_alloced_all, ogmem_nanify
      integer numheap
      integer numcount
      common/gmemdata/iqoff,iq_heap(mxheap),
     &     itag_heap(mxheap), igmem_count, numheap, 
     &     igmem_priority(mxheap),   igamem_priority(mxheap),
     &     igmem_size(mxheap),       igamem_size(mxheap),
     &     igmem_totsize(0:mxcount), igamem_totsize,
     &     igmem_maxsize(0:mxcount), igamem_maxsize, 
     &     igmem_print, numcount,
     &     ogmem_debug, ogmem_alloced_all, ogmem_nanify

      character*16 zgmem_varid,zgamem_arrnam
      common/gmemtext/zgmem_varid(mxheap),zgamem_arrnam(mxheap)
 
      integer mxresreg, mxressect, iadr, nres, nresreg,isize
      parameter(mxresreg=1)
      parameter(mxressect=40)
      common/gmemresv/iadr(0:mxressect,mxresreg),
     &     nres(mxresreg),isize(mxresreg),nresreg
      integer IGMEM_QUIET
      integer IGMEM_NORMAL
      integer IGMEM_DEBUG
      parameter (IGMEM_QUIET =-12)
      parameter (IGMEM_NORMAL=-11)
      parameter (IGMEM_DEBUG =-10)
c     8note* tractlt noo included in vbdens.m mkd1vb
      integer iprint,nbasis,ncol,nsa,ncurt,lenact,lenbas,nactiv
      integer num3,iblkq,ionsec,nbsort,isecdu,n2int,setintpos
      logical incsort,oprs
      integer nnsort,isort,ksort,maxsort
      real*8 scri
c
      common /tractlt/ scri,iprint,nbasis,ncol,nsa,ncurt,lenact,lenbas,
     &                 num3,iblkq,ionsec,nbsort,isecdu,n2int,nactiv,
     &                 incsort,nnsort,isort,ksort,maxsort,oprs
c
c...  incsort : indicates if incore sorting is on
c...  nnsort : # sortbuffers
c...  isort : sortbuffer done
c...  ksort : adres of sortbuffers in qq (igmem_alloc)
c...  maxsort : maximum sort buffer used
c
c...  frozen core definition for TURTLE
c
      real*8 dxyz
      integer ncore,mapcie,mapiee,ncore_i,nsplice
      parameter (nsplice=512)
      common/splice/dxyz(3),ncore,ncore_i,
     1              mapcie(nsplice),mapiee(nsplice)
c
      
      integer ird, iwr, ipu, main, ibl2m, idaf, ibl3d, num8
      integer lenbl, ntapes, iblkci, numlib, iblkpl
      integer numdis, ibldis, irdcpy
      common /iofile/ ird,iwr,ipu,main,ibl2m,idaf,ibl3d,num8,
     +                lenbl(4),ntapes(5),iblkci,numlib,iblkpl,
     +                numdis,ibldis,irdcpy
c
c
c     core :  frozen core energy for the whole job
c     core_i :  frozen core energy  in makeps0 
c
      real*8 core,pot_nuc,core_i
      integer nfil,nofil,iblvb,lblvb,iblhs
      common /ffile/ core,core_i,pot_nuc,nfil,
     1               nofil(20),iblvb(20),lblvb(20),iblhs
      logical ohmat
      integer nstruc,iprinv,ihfile,ihbloc,ipbloc,idumhs
      common/hsinfo/ohmat(3),nstruc,iprinv,ihfile,ihbloc,ipbloc,idumhs
c
c... common to  sitmch from 8 to 16 bits orb ital packing
c...  n8_16 is # bits used for packing; i.e. 8 (255 orbitals) 16 (32000)
c...  n16_32 is the correspnding one for 16 vs 32 bits packin 
c...  GAMESS  parameters cannot be used as VB may have more MO's than AO's
c...  parameters are set in vbcrestr
c
      integer n8_16,n16_32,n340
      common/c8_16vb/n8_16,n16_32,n340
c
c...  find the orthogonality classes in the orbitals
c
      kiort = igmem_alloc_inf(ncore+nsa,'vbci.m','buildnat','kiort',
     &                        IGMEM_DEBUG)
      ks = igmem_alloc_inf(lenact,'vbci.m','buildnat','ks',IGMEM_DEBUG)
c...  skip h
      call rdedx(qq(ks),lenact,iblhs,num3)
      call reads(qq(ks),lenact,num3)
      call sym1(qq(ks),nsa,qq(kiort),northo,'print')
      call gmem_free_inf(ks,'vbci.m','buildnat','ks')

c
c...  read information from datatape
c
c
c.....retrieve size information for memory partitioning
c
      call wr15vb(nelec,nalfa,nconf,ndets,nstruc,maxdet,maxstr,
     &            nwpack,ncoeff,imax,ndum11,ndoub,'read')
c
      kidet = igmem_alloc_inf(nwpi(nelec*ndets),'vbci.m','buildnat',
     &                        'kidet',IGMEM_DEBUG)
      kjdet = igmem_alloc_inf(nwpi(nelec*ndets),'vbci.m','buildnat',
     &                        'kjdet',IGMEM_DEBUG)
      kdete = igmem_alloc_inf(2*ndets,'vbci.m','buildnat','kdete',
     &                        IGMEM_DEBUG)
c
      kidps = igmem_alloc_inf(ndets*nstruc,'vbci.m','buildnat','kidps',
     &                        IGMEM_DEBUG)
      kndet = igmem_alloc_inf(nstruc,'vbci.m','buildnat','kndet',
     &                        IGMEM_DEBUG)
      kcoef = igmem_alloc_inf(ncoeff,'vbci.m','buildnat','kcoef',
     &                        IGMEM_DEBUG)
      kpacd = igmem_alloc_inf(nwpack,'vbci.m','buildnat','kpacd',
     &                         IGMEM_DEBUG)
      kigro = igmem_alloc_inf(5*nstruc,'vbci.m','buildnat','kigro',
     &                        IGMEM_DEBUG)
      kmaxmem = igmem_max_memory()
      kmemscr = kmaxmem/10
      kscr  = igmem_alloc_inf(kmemscr,'vbci.m','buildnat','kscr',
     &                        IGMEM_DEBUG)

      call inform(nelec,nalfa,ndets,nstruc,qq(kpacd),qq(kndet),
     &            qq(kidps),qq(kcoef),ncoeff,qq(kigro),ngroup,
     &            qq(kiort),northo,qq(kscr),kmemscr)

      call gmem_free_inf(kscr,'vbci.m','buildnat','kscr')
c
c...  generate psi0 on determinant basis and unpack determinants
c


      call psi0det(ci,nstruc,qq(kigro),ngroup,qq(kcoef),ncoeff,
     &             ndets,qq(kndet),qq(kidps),qq(kidet),qq(kjdet),
     &             qq(kpacd),nelec,nalfa,ndettot,qq(kdete),kmemscr,
     &             'save')
c
c...  unpack information of mo's per determinant
c
      call unpack(qq(kpacd),n8_16,qq(kidet),ndettot*nelec)
c
      call gmem_free_inf(kigro,'vbci.m','buildnat','kigro')
      call gmem_free_inf(kpacd,'vbci.m','buildnat','kpacd')
      call gmem_free_inf(kcoef,'vbci.m','buildnat','kcoef')
      call gmem_free_inf(kndet,'vbci.m','buildnat','kndet')
      call gmem_free_inf(kidps,'vbci.m','buildnat','kidps')
c
c...  generate natural orbitals
c
      kdmo = igmem_alloc_inf(nsa*(nsa+1)/2,'vbci.m','buildnat',
     &                       'kdmo',IGMEM_DEBUG)
      kdao = igmem_alloc_inf(nbasis*(nbasis+1)/2,'vbci.m','buildnat',
     &                       'kdao',IGMEM_DEBUG)
      ksao = igmem_alloc_inf(nbasis*(nbasis+1)/2,'vbci.m','buildnat',
     &                       'ksao',IGMEM_DEBUG)
      kvnat = igmem_alloc_inf(nbasis*nbasis,'vbci.m','buildnat',
     &                        'kvnat',IGMEM_DEBUG)
      kipos = igmem_alloc_inf(nelec*nelec+2*(nelec*(nelec+1)/2)**2,
     &                        'vbci.m','buildnat','kipos',IGMEM_DEBUG)
c
      kjjdet = igmem_alloc_inf(nwpi(nelec),'vbci.m','buildnat',
     &                        'kjjdet',IGMEM_DEBUG)
      khh = igmem_alloc_inf(lenact,'vbci.m','buildnat','khh',
     &                      IGMEM_DEBUG)
      kss = igmem_alloc_inf(lenact,'vbci.m','buildnat','kss',
     &                      IGMEM_DEBUG)
      kgg = igmem_alloc_inf(lenact*(lenact+1)/2+1,'vbci.m','buildnat',
     &                      'kgg',IGMEM_DEBUG)
c
      call rdedx(qq(khh),lenact,iblhs,num3)
      call reads(qq(kss),lenact,num3)
      call getin2(qq(kgg))
c
      call natorbt(qq(kdete),ndettot,qq(kidet),qq(kjdet),
     &     qq(kjjdet),nelec,nalfa,qq(kdao),qq(kdmo),qq(kipos),
     &     qq(kvnat),qq(kiort),northo,nbasis,nsa,ncore,
     &     qq(khh),qq(kss),qq(kgg+1),qq)
c

      call gmem_free_inf(kgg,'vbci.m','buildnat','kgg')
      call gmem_free_inf(kss,'vbci.m','buildnat','kss')
      call gmem_free_inf(khh,'vbci.m','buildnat','khh')
      call gmem_free_inf(kjjdet,'vbci.m','buildnat','kjjdet')
      call gmem_free_inf(kipos,'vbci.m','buildnat','kipos')
      call gmem_free_inf(kvnat,'vbci.m','buildnat','kvnat')
      call gmem_free_inf(ksao,'vbci.m','buildnat','ksao')
      call gmem_free_inf(kdao,'vbci.m','buildnat','kdao')
      call gmem_free_inf(kdmo,'vbci.m','buildnat','kdmo')
c
      call gmem_free_inf(kdete,'vbci.m','buildnat','kdete')
      call gmem_free_inf(kjdet,'vbci.m','buildnat','kjdet')
      call gmem_free_inf(kidet,'vbci.m','buildnat','kidet')
      call gmem_free_inf(kiort,'vbci.m','buildnat','kiort')
c...  that's it
      return
      end

      subroutine vberr(string)
c
      implicit real*8 (a-h,o-z), integer (i-n)
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
c...     control parasllel working of vb (parallel directive)
      integer sync,check
      common/vbparr/sync,check
c
      character*(*) string
c
c...  vb error routine
c
c
      me = 0
c
      ll = len(string)
      write(iwr,1) me,string(1:ll)
1     format(///,1x,79('*'),/,1x,i4,':',1x,a,/,1x,79('*'),/
     *       /,10x,' no problem is so big or so complicated ...',/
     *        ,10x,' ........... that it can''t be run away from',/      
     *        ,10x,'                               (Schulz)     ',/
     *       /,1x,79('*'))
c
      if ( sync .eq. 0 ) then
        call caserr(" VB error ")
      else if ( sync .eq. 1 ) then
        call caserr2(" VB error ")
      end if
c
      return
      end

      subroutine cjacdav(q,lword,potnuc,energy)
c
c     routine to call jacdav
c     single vector for now
c     Xinyi Xian (june, 2001)
c
      implicit real*8 (a-h,o-z), integer (i-n)
c
      common/davcom/eshift,thresj,thresd,maxcyc,maxdav,maxsel,iprind,
     *              alter,firsts
      logical alter,firsts
      logical ohmat
      integer nstruc,iprinv,ihfile,ihbloc,ipbloc,idumhs
      common/hsinfo/ohmat(3),nstruc,iprinv,ihfile,ihbloc,ipbloc,idumhs
c
      
      integer ird, iwr, ipu, main, ibl2m, idaf, ibl3d, num8
      integer lenbl, ntapes, iblkci, numlib, iblkpl
      integer numdis, ibldis, irdcpy
      common /iofile/ ird,iwr,ipu,main,ibl2m,idaf,ibl3d,num8,
     +                lenbl(4),ntapes(5),iblkci,numlib,iblkpl,
     +                numdis,ibldis,irdcpy
c
c
      dimension q(lword)
c
      ku = 1
      kHu = ku + nstruc
      kSu = kHu + nstruc
      kHd = kSu + nstruc
      kSd = kHd + nstruc
      kt = kSd + nstruc
      kscr1 = kt + nstruc
      kscr2 = kscr1 + nstruc
      kscr3 = kscr2 + nstruc
      kscr4 = kscr3 + nstruc
      kscr5 = kscr4 + 3*nstruc
      kv = kscr5 + 3*nstruc
      kHv = kv + nstruc*maxdav
      kSv = kHv + nstruc*maxdav
      kHtri = kSv + nstruc*maxdav
      kAtri = kHtri + maxdav*(maxdav+1)/2
      kew = kAtri + maxdav*(maxdav+1)/2
      kev = kew + maxdav
      ksi = kev + maxdav*maxdav
      kis = ksi + 2*maxdav
      kend = kis + maxdav
      if (kend.gt.lword) call caserr('core overflow in cjacdav')
c  
      call jacdav(q(ku),q(kHu),q(kSu),q(kHd),q(kSd),q(kt),
     1            q(kscr1),q(kscr2),q(kscr3),q(kscr4),q(kscr5),
     2            q(kv),q(kHv),q(kSv),
     2            q(kHtri),q(kAtri),q(kew),q(kev),q(ksi),q(kis),
     3            nstruc,maxdav,thresd,thresd/10.0d0,energy,
     4            mspace,ihfile,ihbloc,iprind,potnuc,iwr)
c
      end


ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc



       subroutine jacdav(u,Hu,Su,Hd,Sd,t,scr1,scr2,scr3,scr4,scr5,
     +                   v,Hv,Sv,Htri,Atri,ew,ev,si,is,
     +                   n,mmax,epsilon,thresh,energy,mspace,
     +                   ihfile,ihbloc,iprind,potnuc,iwr)
c
c     Reference :
c     G.Sleijpen and H.van der Vorst.Jacobi-Davidson Methods 
c     (Section 5.6). 
c     In Z.Bai,J.Demmel,J.Dongarra,A.Ruhe, and H.van der Vorst, editors,
c     Templates for the Solution of Algebraic Eigenvalue Problems:
c     A Practical Guide, page 123-127, SIAM, Philadelphia, 2000.
c     
c      Initialization phase
c      To apply this algorithm you need to specify a tolerance epsilon, 
c      a target value tau, and a number Kmax that specifies how many 
c      eigenpairs  near tau should be computed. In our case Kmax=1.
c      A starting vector v0 is computed.
       
       implicit none

       integer n,mmax
       real*8 Hd(n),Sd(n),t(n)
c       Hd=diagonal H-matrix; Sd=diagonal S-matrix; t=update vector
       real*8 u(n),Hu(n),Su(n)
c        u= approximate eigenvector; Hu= H-matrix*u; Su= S-matrix*u
c        computed as: Su= Sv*ev and Hu= Hv*ev
       real*8 v(n,mmax),Hv(n,mmax),Sv(n,mmax)
c       v= expansion-space; Hv and Sv= H or S times the expansion-space
       real*8 Atri(*),Htri(*),ev(mmax,mmax),ew(mmax)
c        The interaction-matrix is M=v(i)*H*v(j)
c        Atri= lower diagonal of M deformed in the "jacodiag"
c        Htri= Lower diagonal of M;
c        ev= eigenvector of Hlittle; ew= eigenvalue of Hlittle;

       real*8 tau,epsilon,energy,rr,thresh,aa,ddot,dsqrt,conv,dum
       Integer mspace,i,j,k,m,l,ncyc,is(*),si(*),itri
       Integer ihfile,ihbloc,iprind,iwr
c       (is(n),si(2*n) used in jacobi )
       real*8 scr1(n),scr2(n),scr3(n),scr4(n,3),scr5(n,3),potnuc
c
       itri(i,j) = max(i,j)*(max(i,j)-1)/2 + min(i,j)
c      puts the elements of the matrix in row
c
       mspace  = 0
c
c      dimension of the searchspace is mspace
c
c-----initialisation---------------------------------------------------------
c
       call getd(Hd,Sd,n,ihfile,ihbloc)
c
c....  start-vector = unit vector corresponding to lowest diagonal
c
       tau = 1.0d99
       do i=1,n
        if (Hd(i)/Sd(i).lt.tau) then
          m = i
          tau = Hd(i)/Sd(i)
        end if
       end do
       if (tau.gt.1.0d10) stop " tau too big"
        do i=1,n
         t(i) = 0.0d0
        end do
       t(m) = 1.0d0
c
100    mspace= mspace + 1
c
c       The matrices v and Hv and Sv are 1 enlarged
c
c      orthogonalise t (over S) onto (previous) vectors v by means
c      of Gram-Schmidt.
c      If mspace=1 this is an empty loop
c
       do i=1,mspace-1
          aa= ddot(n,t,1,Sv(1,i),1)
          call daxpy(n,-aa,v(1,i),1,t,1)
       end do
       call hsc(dum,Sv(1,mspace),t,n,ihfile,ihbloc,'s')
       aa = ddot(n,t,1,Sv(1,mspace),1)
       aa = 1.0d0/dsqrt(aa)
c       line (7) from page 124
       call dscal(n,aa,t,1)      
       call dscal(n,aa,Sv(1,mspace),1)
       call dcopy(n,t,1,v(1,mspace),1)
       call hsc(Hv(1,mspace),dum,v(1,mspace),
     +          n,ihfile,ihbloc,'h')
              
c      matrix vector multiply Hv(m)=H*v(m) and Sv(m)=S*v(m) (done above)
c       V denotes the matrix with the current basis vectors v(i) for the
c       search subspace as its columns; likewise Hv and Sv
c      
c       Projection of the H-matrix onto the v-space
c       the matrix is then M=v(i)*H*v(j)
c       add a row to the matrix M (the v's are S-orthonormal)
c
         do l=1,mspace
           Htri(itri(l,mspace)) = ddot(n,v(1,l),1,Hv(1,mspace),1)
         end do
c
c        Copy of Htri into Atri, because the subroutine "jacodiag"
c        changes Htri, but the unchanged version will be needed later
c
         call dcopy(itri(mspace,mspace),Htri,1,Atri,1)
         call jacodiag(Atri,mspace,ev,ew,is,si)
c
c        take eigenvectors with lowest eigenvalue (ew)
c        jacodiag sorts all eigenvectors with the smallest one first
c        so by taking the first one, we have the ev with the lowest ew
         
       energy = ew(1)
c
c      now generate u = approx. eigenvector
c                   Hu and Su
c
       call jd_ev(v(1,1),ev,u,n,mspace)
       call jd_ev(Sv(1,1),ev,Su,n,mspace)
       call jd_ev(Hv(1,1),ev,Hu,n,mspace)
       
c
c
c       residue     r = Hu-energy*Su 
c                     = Hu-ew*Su =(H-energyS)u
c
c      generate r in Hu

       call daxpy(n,-energy,Su,1,Hu,1)
       
c          Stopping criterion:
c          If sqrt(r*r)<eps, the aproximated eigenvector u is ok
c
           rr=ddot(n,Hu,1,Hu,1)
           if (iprind.ge.10) write(iwr,600) mspace,dsqrt(rr),
     1                                      energy+potnuc
600        format(' davidson cycle ',i5,' residue ',e10.3,
     1            ' energy ',f15.8) 
           
           If (dsqrt(rr).lt.epsilon) then
              energy = energy + potnuc
              if (iprind.ge.5) write(iwr,601) mspace,dsqrt(rr)
601           format(' davidson converged in',i6,
     1               ' cycles with residue',e12.5) 
              return
           end if         
           
c          The dimension of the searchspace shouldn't go over mmax
           
           if (m.gt.mmax) then 
            stop
c         
          end if 

c-----------------generation of new update vector t-----------------
c
      call jd_iterate(u,Su,Hu,Hd,Sd,ew,thresh,scr1,scr2,scr3,scr4,scr5,
     +             ncyc,t,n,conv,ihbloc,ihfile,iprind,iwr)
c
c           the update vector t is not really orthogonal to Su
c           because of the preconditioner, which changes the 
c           vector t, so to make sure t is Su-orthogonal, we
c           orthogonalize it again
c
            aa = ddot(n,t,1,Su,1)
c            orthogonalize t and Su
            do i=1,n
               t(i)=t(i)-u(i)*aa
            end do
           
          go to 100
c
       End


ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


      subroutine jd_iterate(u,Su,Hu,Hd,Sd,ew,thresh,
     #                      d,At,Bt,pt,pAt,
     #                      ncyc,t,n,conv,ihbloc,ihfile,iprind,iwr)
c
c         Solve t from C*t=-r with P*t=0 and p=Su so Su*t=0  
c         and C=(I-p*(u*))(H-ew*S)(I-u*(p*))
c         so C*t=-r and r=uA
c            (C+D-D)t=-r with D=diagonal matrix of C
c            (C+D)t-Dt=-r          or  (C-D)t+Dt=-r  
c                  -Dt=-r-(C+D)t           Dt=-r-(C-D)t
c                  -Dt=-r-Ct-Dt            Dt=-r-Ct+Dt
c                   Dt=r+Ct+Dt             Dt=-(r+Ct)+Dt
c                   t=1/d*(r+Ct)+t          t=-1/d*(r+Ct)+t
c         this program is about solving the iterative
c         equation: t'=d*(r+Ct)+t until the res,
c          wich is r+Ct,is below thresh.
c
c         Ct = At -ew.Bt - p.u.At + ew.p.u.Bt
c
c         - input -
c         u : current approx. eigenvector(input)
c         Su : S*u   (S is the metric)  (input)
c         r : is residue of eigenvalue problem
c         Hd : diagonal  of H-matrix
c         Sd : diagonal  of S-matrix
c         ew : current approx. eigenvalue  (input)
c         required is a routine hsc that produces the At,Bt  products
c         d is the  diagonal of C
c         thresh: limit for the res
c
c         - output -
c         t= the best vector of the correction equation
c         ncyc-> number iterations it took (max 50)
c         conv-> convergence

          implicit none
c
          integer n
c
          real*8 u(n),Su(n),Hu(n),Hd(n),Sd(n),thresh,ew
          Integer ncyc,maxcyc
          real*8 t(n),conv,shift,dnormm,scale,convp
c
          real*8 d(n),At(n),Bt(n),pt(n,3),pAt(n,3)
          real*8 uAt,uBt,uBAt,ddot,eu1,eu2,aa

          Integer m,k,i,j,ihbloc,ihfile,iprind,iwr
          real*8 limit
          parameter (limit=1.0d+2,maxcyc=20,scale=0.3d0)
c        
          ncyc=0
          shift = 0.0d0
          dnormm = 1000.0d0
c         call dinicg
c
c----------- initialisation ------------------------------------------
c         
c         use diagonal of (I-p*(u*))(H-ew.S)(I-u*(p*)) as preconditioner
c         the diagonal is:
c 
c         d(i)=Hd(i)-2Au(i)p(i)-ew.Sd(i)+2ew.Su(i)p(i)+eu1u(i)p(i)
c              -ew*eu2u(i)p(i)
c         to spare a vector we call Au =At and Bu =Bt

c---------------------------------------------------------------------

          call hsc(At,Bt,u,n,ihfile,ihbloc,'b')
          eu1= ddot(n,u,1,At,1)         
          eu2= ddot(n,u,1,Bt,1)                  
          
          do i=1,n
          d(i)=Hd(i)+2.0d0*ew*Bt(i)*Su(i)-2.0d0*At(i)*Su(i)-ew*Sd(i)+
     +            Su(i)*Su(i)*eu1-ew*eu2*Su(i)*Su(i)
c                 d(i)=1.0d0/d(i)
c         if (dabs(d(i)).gt.limit) d(i)= dsign(limit,d(i))
          if (dabs(d(i)).lt.1.0d0/limit) d(i)= dsign(1.0d0/limit,d(i))
          end do
          
c         t'=-1/d*(r+Ct)+t
c         start at t=o
c         so the first t'=-r/d   

          do i=1,n
             t(i)=-Hu(i)/(d(i)+shift)
          end do         
                              
*         begin iterations
*
10        ncyc= ncyc +1
*
*         generate Ct = At -o.Bt - p.u.At + o.p.u.Bt
*                      = At -o.Bt + p.(o.u.Bt - u.At)
*         result in At
*
             call hsc(At,Bt,t,n,ihfile,ihbloc,'b')
*         Bt => ew.Bt
             call dscal(n,ew,Bt,1)
             uAt= ddot(n,u,1,At,1)         
             uBt= ddot(n,u,1,Bt,1)         
             uBAt = uBt - uAt
c            
             call daxpy(n,-1.0d0,Bt,1,At,1)
             call daxpy(n,uBAt,Su,1,At,1)
*
*         Ct generated in At
*
*         r+Ct in At
c
             call daxpy(n,1.0d0,Hu,1,At,1)
c
c...   Calling DIIS to prop up residue and solution vector.
c

c       convp = dsqrt(ddot(n,At,1,At,1))
c      call diiscg(t,At,pt,pAt,n)
*
*         r+Ct in At ; check convergence
*
             conv = ddot(n,At,1,At,1) 
             conv = dsqrt(conv)
c             print *, ' rewco convp,conv ',convp,conv
c
c...   Dynamic shifting
c
             dnormm = dmin1(dnormm,conv)
             if (conv.gt.dnormm.and.dnormm.gt.0.0d0) then
                shift = shift + scale*(conv/dnormm - 1.0d0)
             endif
c
             if(iprind.ge.15) write(iwr,603) ncyc,conv,shift      
603          format(' JD-iteration',i6,' residue',e10.2,' shift',e10.2)
             if (conv.lt.thresh) return
             if (ncyc.gt.maxcyc) then
                if (iprind.ge.9) write(iwr,*) ' ncyc gt ',maxcyc,conv
                return
             end if
*          
*         update t : t' = t -(r+Ct)/D 
*
              do i=1,n
               At(i) = At(i)/(d(i)+shift)
              end do
c       (r+Ct)/D => At
             call daxpy(n,-1.0d0,At,1,t,1)
c             
             if(iprind.ge.100) write(iwr,602)t            
602          format('the update t is',/,(1x,10f13.5))
c         
c 

          go to 10
*        
          
          end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc



      subroutine jd_ev(v,s,u,n,ns)
c
c      calculates current approximation of eigenvector
c      or matrix * eigenvector
c
c      v : expansion space
c      s : eigenvector in terms of expansion space
c      u result
c   
       implicit none
c
       integer i,ns,n
       real*8 v(n,ns),s(ns),u(n)
c
c      the starting-vector u is made
c         u = s(1)v(1,i) 
c
       call dcopy(n,v(1,1),1,u,1)
       call dscal(n,s(1),u,1)
c
c       u = u + s(i)v(1,i)
c            
       do i=2,ns
          call daxpy(n,s(i),v(1,i),1,u,1)
       end do     
     
       return
       end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine jacodiag(matrix,n,evecs,evals,is,si)
c
      implicit real*8  (a-h,o-z) , integer   (i-n)
      real*8 matrix
c....
c.... diagonalises triangular matrix with dimension n eigenvectors
c.... are returned in evecs, eigenvalues in evals, is is a scratch
c.... array of dimension n+1, si scratch array 2*n
c
      dimension matrix(*),evecs(*),evals(*),is(*),si(*)
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
c...            symcri = critor*100 remcri = 1.0d-10 crihyb = 1.0d-14
c...            cribri = 1.0d-8 crinor = critor*100  cribub = 1.0d-10
c...            criign = 1.0d-44 cripop 1.0d-10
      real*8 critor,crilow,cridep,symcri,remcri,crihyb,cribri,
     1     crinor,cribub,criign,cripop
      common/vbcri/ critor,crilow,cridep,symcri,remcri,crihyb,cribri,
     1              crinor,cribub,criign,cripop
c
      do i=1,n+1
         is(i)=i*(i-1)/2
      end do
      call jd_jacobi(matrix,is,n,evecs,n,evals,0,2,crilow,si)
      return
      end
      subroutine jd_jacobi(h,iky,nbasis,v,nrow,e,init,iorder,small,y)
c
      implicit real*8  (a-h,o-z) , integer   (i-n)
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
      if (init.ne.1) then
         call vclr(v,1,nrow*nbasis)
         do 10 i=1,nbasis
            v(i,i) = 1.0d0
10       continue
      end if
      if (nbasis.eq.1) then
         e(1)   = h(1)
         return
      end if
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
         return
      else if (iorder.eq.1) then
         call gather(nbasis,e,h,iky(2))
         return
      else if (iorder.eq.2) then
         call gather(nbasis,e,h,iky(2))
         call jd_order(e,v,nbasis,nrow,y,y(nbasis+1), 1)
         call scatter(nbasis,h,iky(2),e)
      else if (iorder.eq.3) then
         call gather(nbasis,e,h,iky(2))
         call jd_order(e,v,nbasis,nrow,y,y(nbasis+1),-1)
         call scatter(nbasis,h,iky(2),e)
      end if
      return
      end

      subroutine igather(n,r,a,map)
      implicit integer (a-z)
      dimension a(*),r(n),map(n)
c
      do 10 loop=1,n
   10 r(loop) = a(map(loop))
c
      return
      end

      subroutine jd_order(e,v,n,nrow,y,imap,ipar)
c
      implicit real*8  (a-h,o-z) , integer   (i-n)
c
      dimension e(*),v(nrow,n),y(*),imap(*)
      do 10 i=1,n
         imap(i) = i
10    continue
      call bubble(e,n,ipar,imap)
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
