c
c  $Author: jvl $
c  $Date: 2010-12-06 15:03:09 +0100 (Mon, 06 Dec 2010) $
c  $Locker:  $
c  $Revision: 6213 $
c  $Source: /c/qcg/cvs/psh/GAMESS-UK/vb/vbutil.m,v $
c  $State: Exp $
c  $Log: not supported by cvs2svn $
c  Revision 1.31  2007/09/25 10:35:51  jvl
c  Fixed pg_sndrcv (call to sendrecv which calls MPI_Sendrecv) because of
c  integer*8/integer*4 type mismatch (for big endians), which lead to
c  checksum error.
c  /marcin
c
c  Revision 1.30  2007/07/24 20:27:51  psh
c  replaced 2 holleriths which crash gfortran compiler
c
c  Revision 1.29  2007/07/06 18:45:44  jvl
c  Minor bug fix in CRESTR with STRUC NORM.
c  Minor update of VBMO directive.
c  Addition of new stripblanks subroutine to strip all the front and trailing
c  blank spaces from string, while keeping all the within.
c  /marcin
c
c  Revision 1.28  2007/07/04 12:37:56  jvl
c  Introducing new VB switch - VBMO, first step to separate output file for molden.
c  Now, upon request, separate file with VB orbitals (with proper format) can be created,
c  called out.vbmo.<title_from_VBMO_switch>. Check vbin.m comment/GAMESS manual for more.
c  Major bugfixes of dimensions in makepsi0 subroutine:
c   - serial + CORE XX
c   - parallel + bigger basis set than 6-31g + CORE switch
c  Minor changes in igmem_alloc_inf names.
c  /marcin
c
c  Revision 1.27  2007/03/20 14:49:31  jvl
c  Pretty major overhaul of the VB code
c  say 50% is now dynamic using igmem_alloc and, oh wonder,
c  all the examples are still checking out OK, including a new resonating
c  super hybrid. (will follow)
c
c  Revision 1.26  2006/04/18 16:35:34  jvl
c  - added basis option to allow Vb to be run in anothr basis.
c  - cleaned up various bits of VB, like core partitioning
c  - adapted get1e
c
c  Revision 1.25  2006/01/13 17:56:48  jmht
c  Merging in changes from the release branch
c
c  Revision 1.24  2005/09/23 14:34:56  jmht
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
c  Revision 1.23.2.1  2005/07/19 06:52:20  wab
c
c  min0/max0 changes (see m4 checkin)
c
c  Revision 1.23  2005/05/11 22:37:33  jvl
c  made paralleol vb - I8 completer + few checks and better errormessages
c  sendrecv introduced to avoid race conditions (added to ga-stuff really)
c
c  Revision 1.22  2005/04/22 15:07:55  jvl
c  added mpi_sndrcv to get vb 4-index going
c
c  Revision 1.21  2005/04/22 11:07:56  jvl
c  All vb criteria are now in common vbcri and may be set
c  in the input. Also their defaults are now consistent
c
c  Revision 1.20  2005/03/25 23:16:58  jvl
c  changes to allow I8 VB; changing individual integer labels to packed reals
c
c  Revision 1.19  2005/02/05 18:04:26  wab
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
c  Revision 1.18  2004/03/24 14:57:46  jvl
c  fixed a bug in print-out of branching diagrams and leading terms (JJE).
c
c  Revision 1.17  2003/11/30 14:50:41  jvl
c  - made makefile to make itanium with efc
c    also added -ftz to flush to 0.0 and -cm to restricyt stupid warnings
c  - replaced flush by flushout for unit 6 for altix (flush comes
c    from a library and is something else ???)
c  - fixed a few complaints from efc compiler
c  - move c2030_b and c2030_c (reluctantly) to vb
c    left the files also in chap2 but commented call and check
c  Now single processor I4 version is fully (??) OK on ALTIX
c
c  Revision 1.16  2003/10/08 17:31:14  jvl
c  several bugfixes (JJE 8/10/2003).
c
c  Revision 1.15  2003/08/27 12:55:32  jvl
c  fixed usage of non-standard fortran with parameters (JJE 27/8/2003).
c
c  Revision 1.14  2003/08/27 10:40:34  jvl
c  bugfix in prweig (JJE 27/8/2003)
c
c  Revision 1.13  2003/08/26 13:36:00  jvl
c  Weights are printed together with their branching or rumer diagram.
c  (JJE 26/8/2003)
c
c  Revision 1.12  2003/06/10 16:43:03  jvl
c  Minor bugfixes (JJE 10/6/2003)
c
c  Revision 1.11  2003/02/18 17:17:18  jvl
c  A lot of changes, in general:
c  - In vbci the call to bsmat for natural orbitals is no longer necessary
c  - The suboptions pert and fock for optimise in scf are implemented
c  - The subroutine vbscf is completely restructured and is more logical
c  - Addition of the makepsi0, which uses a small extra 2-electron,
c    transformation
c    JJE - 18 feb 2003
c
c  Revision 1.10  2003/01/19 16:23:53  jvl
c  attempt to correct a comment (in the previous cvs checking) that was interpreted
c  by m4 .....; hope a fixed it ......
c
c  Revision 1.9  2003/01/19 15:52:35  jvl
c  corrected ifdef parallel 
c  vb now end properly for hcpu
c
c  Revision 1.8  2003/01/19 12:33:07  jvl
c  - check on # active orbitals (an infinite loop could result)
c  - hcpu option to allow primitive restarting of matrix element calculation in vbci
c
c  Revision 1.7  2002/09/05 14:50:32  jvl
c  Changed dimension of igroup from 3 to 5 (has to do with new opti option) (JJE).
c
c  Revision 1.6  2001/10/24 21:45:40  jvl
c  corrected few problems in geometry optimisation with VB and added initial stuff
c  to limit output later..
c
c  Revision 1.5  2001/07/03 14:13:20  jvl
c  fixed a few bugs in parallel code and allowed non-usage of peigss in vb
c  (i.e. peigss specifies if peigas is to be used in vb)
c
c  Revision 1.4  2001/06/12 12:21:21  jvl
c  - Removed several bugs
c  - Improved parallel 4-index transformation including changes in getin2
c  - Added timing for VB routines
c
c  Revision 1.3  2001/02/08 14:34:40  jvl
c  Adapted parallelisation for GAMESS-UK.
c  An MPI and GA version have been created
c
c  Revision 1.2  2000/03/31 14:18:28  jvl
c  Changed vnorm to vnrm
c  Changed get1e to use getmat
c  Changed scopy to icopy
c  Improved on information about titles
c  Fokke Dijkstra
c
c  Revision 1.1  2000/02/16 12:20:48  jvlK
c  adding vb files to gamess repository
c
c Revision 1.13  1998/02/25  13:04:15  fokke
c removed routine flushbuffer (now in util_atmol.F) and inserted common
c block io
c
c Revision 1.13  1998/02/25  13:04:15  fokke
c removed routine flushbuffer (now in util_atmol.F) and inserted common
c block io
c
c Revision 1.12  1998/02/16  15:09:00  fokke
c added mp_kill
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
c Revision 1.5  1996/11/07  17:10:01  joop
c SiGr => SGI + bugfixje
c
c Revision 1.4  1996/11/04  16:03:51  joop
c oipsc,oroot,mypid,iipsc moved to util_atmol.F in lib
c
c Revision 1.3  1996/10/29  14:47:28  joop
c rcs + parallel/labeled matrix elements
c
c
       subroutine abmax(a,n,is,shape,nn,amax,imax)
       implicit REAL (a-h,o-z), integer (i-n)
       dimension a(*)
       character*(*) shape
       itri(i,j)=max(i,j)*(max(i,j)-1)/2+min(i,j)
c
       amax = -1.0d0
       do 10 i=is,n
          if (shape.eq.'tri') then
             ii = itri(i,nn)
          else if (shape.eq.'sq') then
             ii = nn*(i-1)+1
          else
             call caserr('unknown shape in abmax')
          end if
          if (dabs(a(ii)).gt.amax) then
             amax = dabs(a(ii))
             imax = i
          end if
10     continue
c
       return
       end
       subroutine abmaxsel(a,n,is,isel,nsel,shape,nn,amax,imax)
       implicit REAL (a-h,o-z), integer (i-n)
       dimension a(*),isel(nsel)
       character*(*) shape
       itri(i,j)=max(i,j)*(max(i,j)-1)/2+min(i,j)
c
       amax = -1.0d0
       imax = 0 
       do 10 i=is,n
          if (shape.eq.'tri') then
             ii = itri(i,nn)
          else if (shape.eq.'sq') then
             ii = nn*(i-1)+1
          else
             call caserr('unknown shape in abmax')
          end if
          do j=1,nsel
             if (i.eq.isel(j)) go to 5
          end do
          if (dabs(a(ii)).gt.amax) then
             amax = dabs(a(ii))
             imax = i
          end if
    5     continue 
10     continue
c
       return
       end
       subroutine abmin(a,n,amin,imin)
       implicit REAL (a-h,o-z), integer (i-n)
       dimension a(*)
       imin = 1
       amin = dabs(a(1))
       do 10 i=2,n
          if (dabs(a(i)).lt.dabs(amin)) then
             amin = a(i)
             imin = i
          end if
10     continue
       return
       end
      subroutine brtsbh(ibll,ifil)
c
      implicit REAL  (a-h,o-z) , integer   (i-n)
c
c.....
c.....initialise the process of writing a contiguous set of data to
c.....ed(ifile+1),starting at block ibll. this routine is used for
c.....h-matrix elements (brtsbs is used for the s-matrix elements).
c.....
      common/diskh/hof(170),sof(170),cricr(170),no,iblock,ifile,idum,
     +             icof(170),irof(170)
      ifile = ifil
      iblock = ibll
      no = 0
      call search(iblock,ifile)
      return
      end
c
c...  Returns first and last characters, non blank, from string.
c...  Thus, keeps the blanks within string, but only got rid of the
c...  front and trailing ones.
      subroutine stripblanks(fchr,lchr,string)
c
      character string*(*)
      integer i,j,fchr,lchr,stpoint
c
      stpoint = len(string)
      do i=1,stpoint
        if ( string(i:i) .eq. ' ' ) then
        else
          fchr = i
          go to 101
        end if
      end do
      
101   continue

      do j=stpoint,1,-1
        if ( string(j:j) .eq. ' ' ) then
        else 
          lchr = j
          go to 102
        end if
      end do 
      
102   continue

      return
      end
c
      subroutine bubble(a,n,ipar,imap)
c
      implicit REAL  (a-h,o-z) , integer   (i-n)
c
c.....
c.....simple bubble sort. imap keeps track of the permutations
c.....ipar enables one to use this routine for sorting increasingly
c.....(ipar > 0 => first element smallest) or the other way around
c.....
      dimension a(n),imap(n)
      logical ready
10    ready = .true.
      do 20 i=2,n
         if (a(i-1)*ipar.gt.a(i)*ipar) then
            r         = a(i-1)
            a(i-1)    = a(i  )
            a(i  )    = r
            ready     = .false.
            ii        = imap(i-1)
            imap(i-1) = imap(i  )
            imap(i  ) = ii
         end if
20    continue
      if (.not.ready) goto 10
      return
      end
      subroutine bubblevec(a,n,v,nrow,ipar,imap)
c
      implicit REAL  (a-h,o-z) , integer   (i-n)
c
c.....
c.....simple bubble sort. imap keeps track of the permutations
c.....ipar enables one to use this routine for sorting increasingly
c.....(ipar > 0 => first element smallest) or the other way around
c.....when values of e are equal the sort is done on highest vector 
c.....element. when these are equal the element with highes position
c.....comes first (this was done to make the sort reproducable).
      dimension a(n),imap(n),v(nrow,n)
INCLUDE(common/vbcri)
      logical ready
10    ready = .true.
      do 20 i=2,n
         if ((a(i-1)*ipar-a(i)*ipar).gt.cribub) then
            r         = a(i-1)
            a(i-1)    = a(i  )
            a(i  )    = r
            ready     = .false.
            ii        = imap(i-1)
            imap(i-1) = imap(i  )
            imap(i  ) = ii
         else if ((a(i-1)*ipar-a(i)*ipar).gt.-cribub) then
            j1 = imap(i-1)
            j2 = imap(i)
            call maxvec(v(1,j1),n,vmax1,ipos1)
            call maxvec(v(1,j2),n,vmax2,ipos2)
            if ((vmax2-vmax1).gt.cribub) then
              r         = a(i-1)
              a(i-1)    = a(i  )
              a(i  )    = r
              ready     = .false.
              ii        = imap(i-1)
              imap(i-1) = imap(i  )
              imap(i  ) = ii
            else if ((vmax2-vmax1).gt.-cribub) then
               if (ipos1.gt.ipos2) then
                 r         = a(i-1)
                 a(i-1)    = a(i  )
                 a(i  )    = r
                 ready     = .false.
                 ii        = imap(i-1)
                 imap(i-1) = imap(i  )
                 imap(i  ) = ii
               end if
            end if
         end if
20    continue
      if (.not.ready) goto 10
      return
      end
      subroutine ibubdet(ia,n,idir,amap,ipar)
c
      implicit REAL (b-h,o-z), integer (i-n), character (a)
c
c.....
c.....simple bubble sort. amap keeps track of the permutations of spins
c.....idir enables one to use this routine for sorting increasingly
c.....(idir > 0 => first element smallest) or the other way around
c.....
      dimension ia(n),amap(n)
      logical ready
      ipar = 1
10    ready = .true.
      do 20 i=2,n
         if (ia(i-1)*idir.gt.ia(i)*idir) then
            ir         = ia(i-1)
            ia(i-1)    = ia(i  )
            ia(i  )    = ir
            ipar      = -ipar
            ready     = .false.
            aa        = amap(i-1)
            amap(i-1) = amap(i  )
            amap(i  ) = aa
         end if
20    continue
      if (.not.ready) goto 10
      return
      end
      subroutine maxvec(v,n,vmax,ipos)
c
      REAL v(n),vmax
      integer ipos,n
c
      vmax = dabs(v(1))
      ipos = 1
      do i=2,n
        if (vmax.lt.dabs(v(i))) then
          vmax = dabs(v(i))
          ipos = i
        end if
      end do 
      return
      end        
      subroutine writh(iband,igroup,hamil,overl,ic,ir)
c
      implicit REAL  (a-h,o-z) , integer   (i-n)
c
c.....
c.....this routine controls the writing of the h-matrix
c.....ham1 contains the original matrix elements (ordered as groups of
c.....elements). ham2 will contain the same elements, but then reordered
c.....so that, per structure, all matrix elements involving this
c.....structure up until the structure itself appear after one another.
c.....now the lements are written with labels (1996)
c.....
      dimension igroup(5,*),hamil(*),overl(*),ic(*),ir(*)
c.....
c.....originally :
c.....determine the position vector that defines the reordering scheme
c.....now:
c.....determine row and colum indices of nonzero matrix elements
c.....
INCLUDE(common/vbcri)
c
c...  quick and dirty fix // Fokke please consider !!
c
      common/subinf/nwidth
      nwidth = max(nwidth,igroup(2,iband)) 
c
      iold = 1
      inew = 0
      ir0 = 0
      do 5 i=1,iband-1
5     ir0 = ir0 + igroup(2,i)
      ic0 = 0
c
      do 40 k=1,iband
         do 30 irr=1,igroup(2,iband)
            iccend = igroup(2,k)
            if (k.eq.iband) iccend = irr
            do 20 icc=1,iccend
               if (dabs(hamil(iold)).gt.criign) then
                  inew = inew + 1
                  ic(inew) = ic0 + icc
                  ir(inew) = ir0 + irr
                  hamil(inew) = hamil(iold)
                  overl(inew) = overl(iold)
               end if
20          iold = iold + 1
30       continue
         ic0 = ic0 + igroup(2,k)
40    continue
c.....
c.....write to disk
c.....
      call wrtsbh(hamil,overl,ic,ir,inew)
c
      return
      end
      subroutine wrtsbh(h,s,icol,irow,nword)
c
      implicit REAL  (a-h,o-z) , integer   (i-n)
c
c.....
c.... write full records (see brtsbh) this routine is ment for
c.....writing labeled hamiltonian matrix elements
c.....
      dimension h(nword),s(nword),icol(nword),irow(nword)
      common/diskh/hof(170),sof(170),cricr(170),no,iblock,ifile,idum,
     +             icof(170),irof(170)
      integer iprhs,irowhs
      common/hcpu_vb/ iprhs,irowhs
_IF(peigss)
INCLUDE(common/gainfo)
INCLUDE(../m4/common/vcore)
_ENDIF
c
      if (nword.le.0) return
c
      if (iprhs.gt.0) call hspr_long(h,s,icol,irow,nword)
c
_IF(peigss)
      nc = icol(nword)-icol(1)+1
      nr = irow(nword)-irow(1)+1
      ih = igmem_alloc(nc*nr)
      is = igmem_alloc(nc*nr)
      call vclr(Q(ih),1,nc*nr)
      call vclr(Q(is),1,nc*nr)
      do i=1,nword
        ic = icol(i) - icol(1) + 1
        ir = irow(i) - irow(1) + 1
        ipos = ir + (ic-1)*nr -1
        call flushbuffer()
        if (icol(i).ne.irow(i)) then
           Q(ih+ipos) = 2*h(i)
           Q(is+ipos) = 2*s(i)
        else
           Q(ih+ipos) = h(i)
           Q(is+ipos) = s(i)
        endif
      end do
      call pg_acc(iga_h,icol(1),icol(nword),
     &            irow(1),irow(nword),Q(ih),nc,1.0d0)
      call pg_acc(iga_s,icol(1),icol(nword),
     &            irow(1),irow(nword),Q(is),nc,1.0d0)
      call gmem_free(is)
      call gmem_free(ih)
_ELSE
      k = nword
      j = 1
c.....
c.... fill hof
10    nn = min(k,170-no)
      call fmove(h(j),hof(no+1),nn)
      call fmove(s(j),sof(no+1),nn)
      call icopy(nn,icol(j),1,icof(no+1),1)
      call icopy(nn,irow(j),1,irof(no+1),1)
      no = no + nn
      if (no.lt.170) return
c.... buffer full => write out / check positioning
      call search(iblock,ifile)
      call pack(cricr(1),32,icof,170)
      call pack(cricr(86),32,irof,170)
      call put(hof,511,ifile)
      iblock = iblock + 1
      no = 0
      k = k - nn
      j = j + nn
c
      if (k.gt.0) go to 10
_ENDIF
c
      return
      end
      subroutine hspr_long(h,s,icol,irow,nword)
c
      implicit REAL (a-h,o-z)
      dimension h(nword),s(nword),icol(nword),irow(nword)
c
INCLUDE(common/ffile)
INCLUDE(../m4/common/iofile)
c
      logical oroot
      common/junk/ hdum(500),sdum(500)
c
      if (nword.gt.500) call vberr('enlarge routine hspr_long')
      call dcopy(nword,h,1,hdum,1)
      call dcopy(nword,s,1,sdum,1)
_IF(parallel)
      call pg_dgop(253,hdum,nword,'+')
      call pg_dgop(254,sdum,nword,'+')
_ENDIF
      if (oroot()) then
         if (icol(1).eq.1.and.irow(1).eq.1) then
            write(iwr,'(a,e25.16)') ' The CORE contribution is ',core
         end if
         do i=1,nword
            write(iwr,1) irow(i),icol(i),hdum(i)+sdum(i)*core,sdum(i)
1           format(' i ',i3,'  j ',i3,'  h ',e25.16,'  s ',e25.16)
         end do
         call flushbuffer()
      end if
c
      return
      end  
      subroutine corfait(need,navail,string)
c
      implicit REAL (a-h,o-z), integer (i-n)
INCLUDE(../m4/common/iofile)
c
c...  vb core-failure routine
c
      character*(*) string
c
      write(iwr,1) string,need,navail
1     format(//,1x,79('*'),/,' core failure in ',a20,/,
     *       ' we need :',i10,9x,' we have :',i10,/,1x,79('*'))
      call vberr(' core failure ')
c
      return
      end
      subroutine ertsbh
c
      implicit REAL  (a-h,o-z) , integer   (i-n)
c
c.....
c.....close writing of a contiguous array (hamiltonian matrix elements)
c.....
      common/diskh/hof(170),sof(170),cricr(170),no,iblock,ifile,idum,
     +             icof(170),irof(170)
c.... output last (partial) buffer
c.....search must be called first, because overlap matrix elements are
c.....written at the same time
c.....
_IF(peigss)
      return
_ENDIF
      call search(iblock,ifile)
      if (no.ne.0) then
         call pack(cricr(1),32,icof,170)
         call pack(cricr(86),32,irof,170)
         call put(hof,511,ifile)
         iblock = iblock + 1
      end if
      call put(hof,0,ifile)
      iblock = iblock + 1
      no = 0
      return
      end
      subroutine whersbt(ibl,nno)
c
c.... return where we are in wrtsb writing
c.... (i.e. what is next block to be written
c
      implicit REAL  (a-h,o-z) , integer   (i-n)
      common/diskh/hof(170),sof(170),cricr(170),no,iblock,ifile,idum,
     +             icof(170),irof(170)
c
      ibl = iblock
      nno = no
c
      return
      end
      subroutine getd(dh,ds,nstruc,ihfile,ihbloc)
c
      implicit REAL  (a-h,o-z) , integer   (i-n)
c
c...  changed for labeled ints
c
      dimension dh(nstruc),ds(nstruc)
      common/matrixt/hsof(170,2),cricr(170),no,iblock,
     +               icof(170),irof(170)
      character*1 type
_IF(peigss)
INCLUDE(common/gainfo)
_ENDIF
_IF(linux)
      external fget
_ENDIF   
c
      call vclr(dh,1,nstruc)
      call vclr(ds,1,nstruc)
_IF(peigss)
      do i=1,nstruc
         call pg_get(iga_h,i,i,i,i,dh(i),1)
         call pg_get(iga_s,i,i,i,i,ds(i),1)
      end do  
_ELSE
      call search(ihbloc,ihfile)
10    call fget(hsof,nw,ihfile)
      if (nw.eq.0) go to 30
       call unpack(cricr(1),32,icof,170)
       call unpack(cricr(86),32,irof,170)
      do 20 i=1,no
         if (icof(i).eq.irof(i)) then
            if (icof(i).gt.nstruc) call vberr('diag overflow in getd')
            dh(icof(i)) = hsof(i,1)
            ds(icof(i)) = hsof(i,2)
         end if
20    continue
      go to 10
c
30    continue
c
_IF(parallel)
      call pg_dgop(480,dh,nstruc,'+')
      call pg_dgop(4190,ds,nstruc,'+')
_ENDIF
_ENDIF
      end
      subroutine ibubbl(ia,n,ipar)
c
      implicit REAL  (a-h,o-z) , integer   (i-n)
c
c.....
c.....simple integer bubble sort.
c.....ipar enables one to use this routine for sorting increasingly
c.....(ipar > 0 => first element smallest) or the other way around
c.....
      dimension ia(n)
      logical ready
10    ready = .true.
      do 20 i=2,n
         if (ia(i-1)*ipar.gt.ia(i)*ipar) then
            iiii    = ia(i-1)
            ia(i-1) = ia(i  )
            ia(i  ) = iiii
            ready   = .false.
         end if
20    continue
      if (.not.ready) goto 10
      return
      end
      function locati(i,ii,n)
      implicit REAL (a-h,o-z)
      dimension ii(*)
      do 10 j=1,n
         if (ii(j).eq.i) then
            locati = j
            return
         end if
10    continue
      locati = 0
      return
      end
      subroutine makes(s,ir,nr,ic,nc,ipos,super)
c.....
c.....in this routine s is constructed from super using the labels from
c.....ir(ow) and ic(olumn). super is assumed to be a triangular matrix.
c.....
c
      implicit REAL  (a-h,o-z) , integer   (i-n)
c
      dimension s(nr,nc),ir(nr),ic(nc),ipos(nr*nc),super(*)
      common /posit/ iky(3)
      ind(i,j) = iky(max(i,j)) + min(i,j)
      it = 0
         do 20 j=1,nc
            do 10 i=1,nr
               it = it + 1
               ipos(it) = ind(ir(i),ic(j))
10          continue
20       continue
      call gather(it,s,super,ipos)
      return
      end
      subroutine maket(s,v,n,nbasis,sao,scr)
      implicit REAL (a-h,o-z), integer(i-n)
c.....
c.....construct (partial) metric
c.....      
      dimension s(*),v(nbasis,n),sao(*),scr(*)
INCLUDE(common/vbcri)
      it = 0
      do 20 j=1,n
         call cntrc(sao,v(1,j),scr,nbasis,crilow)
         do 10 i=1,j
            it = it + 1
            s(it) = ddot(nbasis,v(1,i),1,scr,1)
10       continue
20    continue
      return
      end
      subroutine matmpl(v,t,vi,n,m)
c
      implicit REAL  (a-h,o-z) , integer   (i-n)
c
c
c     matmpl performs the transformation of a n*m vector-matrix v with
c     a n*n transformation matrix t, both lineairly stored.
c     results are returned in array v.
c     one intermidiate array vi of length n is required
c
      dimension v(1),t(1),vi(1)
      jend = n*m
      do 30 l=1,m
      k = 0
      do 20 i=1,n
      res = 0.0d0
      do 10 j=l,jend,m
      k = k + 1
   10 res = res + t(k)*v(j)
   20 vi(i) = res
      k = 0
      do 30 j=l,jend,m
      k = k + 1
   30 v(j) = vi(k)
c
      return
      end
      subroutine mult11(h,r,scr,nstruc,ndets,q,qdagg)
c
      implicit REAL  (a-h,o-z) , integer   (i-n)
c
c.....
c.....in this routine the matrix-product (q-dagger) * h * q is formed
c.....h is assumed to be a triangular matrix
c.....
      dimension h(*),r(nstruc*nstruc),scr(ndets*nstruc+nstruc*nstruc),
     &          q(ndets*nstruc),qdagg(nstruc*ndets)
c.....
c.....scale the diagonal of h by 0.5
c.....
      do 10 i=1,ndets
         h(i*(i+1)/2) = h(i*(i+1)/2) * 0.5d0
10    continue
      if (ndets/2.lt.nstruc) then
c.....
c.....   choose algorithm 2 (c.f.four index transformation of atmol)
c.....
c.....   zeroize scratch arrays: scr(1 to nstruc*ndets)     = x
c.....                           scr(nstruc*ndets+1 to end) = y
c.....
         call vclr(scr,1,nstruc*ndets + nstruc*nstruc)
c.....
c.....   form x = (q-dagger) * (h-triangle)
c.....
         call dagger(ndets,nstruc,q,ndets,qdagg,nstruc)
         call mxmtturtle(qdagg,nstruc,h,scr,nstruc,nstruc,ndets)
c.....
c.....   form y = x * q
c.....
         call mxmb(scr,1,nstruc,
     &             q,1,ndets,
     &             scr(nstruc*ndets+1),1,nstruc,
     &             nstruc,ndets,nstruc)
c.....
c.....   r(esult-triangle) = y + (y-dagger)
c.....
         call symm1(r,scr(nstruc*ndets+1),nstruc)
c.....
      else
c.....
c.....   choose algorithm 1
c.....
c.....   zeroise scratch space
c.....
         call vclr(scr,1,ndets*nstruc)
c.....
c.....   form h(triangle) * q
c.....
         n = 1
         l = 1
         do 30 i=1,nstruc
            m = 1
            do 20 j=1,ndets
               call daxpy(j,q(n),h(m),1,scr(l),1)
               n = n + 1
               m = m + j
20          continue
         l = l + ndets
30       continue
c.....
c.....   r(triangle) = q(dagger) * scr + scr(dagger) * q
c.....
         m = 1
         n = 1
         do 50 i=1,nstruc
            l = 1
            do 40 j=1,i
               r(m) = ddot(ndets,q(n),1,scr(l),1) +
     &                ddot(ndets,q(l),1,scr(n),1)
               m = m + 1
               l = l + ndets
40          continue
            n = n + ndets
50       continue
      end if
      return
      end
      subroutine mxmtturtle(a,mrowa,b,r,mrowr,ncol,nrow)
c
      implicit REAL  (a-h,o-z) , integer   (i-n)
c
c.....
c.....r = a * b(triangle)
c.....
      dimension a(mrowa,*),b(*),r(mrowr,*)
      m=1
      do 40 i=1,nrow
         do 30 j=1,i
            top = b(m)
            if (top.ne.0.0d0) then
c              do 5 loop=1,ncol
c5             r(loop,i)=r(loop,i)+a(loop,j)*top
               call daxpy(ncol,top,a(1,j),1,r(1,i),1)
            end if
30       m = m + 1
40    continue
      return
      end
       subroutine scalev(a,scal,n)
       implicit REAL (a-h,o-z), integer (i-n)
       dimension a(*)
       do 10 i=1,n
10     a(i) = a(i) * scal
       return
       end

      subroutine tripriacc(a,ndim,isect)
c.....
      implicit REAL (a-h,o-z) , integer (i-n)
INCLUDE(../m4/common/iofile)
c.....
c.....this routine prints (a maximum of 10 columns
c.....on one row) the triangular matrix a
c.....
      dimension a(*)
c.....
      write(iwr,'(/)')
      kend=ndim/4+1
      jstart=1
      it=0                    
      nextra=0
c.....
      do 30 k=1,kend
         n=0
         do 20 j=jstart,ndim
            n=n+1
            if (n.gt.4) then
               iend=4
            else
               iend=n
            end if
            write(isect,10) j,(a(it+i),i=1,iend)
10          format(2x,i3,2x,4f22.16)
            it=it+n+nextra
20       continue
      write(isect,'(/)')
c.....noe=number of elements in a if ndim=10*k
      noe=5*(10*k*k+k)
      nextra=nextra+10
      it=noe+nextra
      jstart=jstart+10
30    continue
c
      return
      end


      subroutine tripri(a,ndim)
c.....
      implicit REAL (a-h,o-z) , integer (i-n)
INCLUDE(../m4/common/iofile)
c.....
c.....this routine prints (a maximum of 10 columns
c.....on one row) the triangular matrix a
c.....
      dimension a(*)
c.....
      write(iwr,'(/)')
      kend=ndim/10+1
      jstart=1
      it=0                    
      nextra=0
c.....
      do 30 k=1,kend
         n=0
         do 20 j=jstart,ndim
            n=n+1
            if (n.gt.10) then
               iend=10
            else
               iend=n
            end if
            write(iwr,10) j,(a(it+i),i=1,iend)
10          format(2x,i3,2x,10f15.7)
            it=it+n+nextra
20       continue
      write(iwr,'(/)')
c.....noe=number of elements in a if ndim=10*k
      noe=5*(10*k*k+k)
      nextra=nextra+10
      it=noe+nextra
      jstart=jstart+10
30    continue
c
      return
      end
      subroutine tripr2(a,ndim)
c
      implicit REAL  (a-h,o-z) , integer   (i-n)
INCLUDE(../m4/common/iofile)
c
c.....
c.....this routine prints the triangular matrix a
c.....
      dimension a(*)
      if (ndim.le.0) write(iwr,20) ndim,a(1)
      it = 0
      do 10 j=1,ndim
         write(iwr,20) j,(a(it+i),i=1,j)
         it = it + j
10    continue
20    format(i6,(t8,10e12.4))
      return
      end

      subroutine prenres(a,ndim,evb)
c
      implicit REAL  (a-h,o-z) , integer   (i-n)
INCLUDE(../m4/common/iofile)
c
c...
c...  print the resonance energy of the molecule
c...
      dimension a(*)
      REAL evb,elow,eres
c
      REAL prlev
      common/leading/prlev
c
      it = 1
      elow = 0.0d0
      do j=1,ndim
         if ( a(it) .lt. elow ) then
           elow = a(it)
         end if
         it = it + j + 1
      end do
      write(iwr,'(a,F14.8)') ' Lowest energy of the structure: ',elow
      eres = (evb - elow)*627.509
      write(iwr,'(a,F8.2,a)') ' The resonance energy: ',eres,' kcal/mol'

      return
      end

      subroutine prdiag(a,ndim,weightgn)
c
      implicit REAL  (a-h,o-z) , integer   (i-n)
INCLUDE(../m4/common/iofile)
c
c...
c...  print the energies of the structures
c...
      dimension a(*)
      REAL weightgn(*)
c
      REAL prlev
      common/leading/prlev
c
      write(iwr,'(/)')
      it = 1
      do j=1,ndim
         if ( weightgn(j) .gt. prlev ) then
           write(iwr,'(3x,I5,2x,F16.10)') 
     &     j,a(it)
         end if
         it = it + j + 1
      end do

      return
      end

      subroutine tripr3(a,ndim)
c
      implicit REAL  (a-h,o-z) , integer   (i-n)
INCLUDE(../m4/common/iofile)
c
c.....
c.....this routine prints the triangular matrix a
c.....
      dimension a(*)
      if (ndim.le.0) write(iwr,20) ndim,a(1)
      it = 0
      write(iwr,'(/)')
      do 10 j=1,ndim
         write(iwr,20) j,(a(it+i),i=1,j)
         it = it + j
10    continue
20    format(3x,i2,(t8,9f12.6))
      return
      end
      subroutine trisqu(a,b,ndim)
c
      implicit REAL  (a-h,o-z) , integer   (i-n)
c
c.....
c.....the triangle a fills the square b
c.....
      dimension a(*),b(ndim,ndim)
      it = 0
      do 20 i=1,ndim
         do 10 j=1,i
            it = it + 1
            b(i,j) = a(it)
            b(j,i) = a(it)
10       continue
20    continue
      return
      end
      subroutine trtrtr(s,p,r,ndim,scr)
c
      implicit REAL (a-h,o-z), integer (i-n)
c
c...  s and p are hermitian matrices in triangular form.
c...  s' * p' * s' is returned (that is the "true" matrix product)
c...  this is not ment to be fast / r may overlap p
c
      dimension s(*),p(*),r(*),scr(ndim,ndim)
      ind(i,j) = max(i,j)*(max(i,j)-1)/2 + min(i,j)
      call vclr(scr,1,ndim*ndim)
      do 30 i=1,ndim
         do 20 j=1,ndim
            do 10 k=1,ndim
               scr(i,j) = scr(i,j) + p(ind(i,k)) * s(ind(k,j))
10          continue
20       continue
30    continue
      call vclr(r,1,ndim*(ndim+1)/2)
      do 60 i=1,ndim
         do 50 j=1,ndim
            do 40 k=1,ndim
               r(ind(i,j)) = r(ind(i,j)) + s(ind(i,k)) * scr(k,j)
40          continue
50       continue
60    continue
      do 80 i=2,ndim
         do 70 j=1,i-1
            r(ind(i,j)) = r(ind(i,j))/2
70       continue
80    continue
      return
      end
      subroutine flushbuffer
c
c...  flush output buffer
c
INCLUDE(../m4/common/iofile)
_IF(sgi)
      call flush(iwr)
_ENDIF
c
      return
      end

_IF(parallel)
      subroutine poe_next()
INCLUDE(../m4/common/iofile)
      common/charsort/sortname(2),scratch,zout
      character*8 zout
      character*44 sortname,scratch
      logical oroot
      character*4 ppp
c
      me = ipg_nodeid()
      if (zout.eq.'output') then
         if (.not.oroot()) open(6,file='/dev/null')
      else if (zout.eq.'ibm') then
c...     ibm may redirect it's output itself
      else 
         if (me.le.9) then
            write(ppp,'(i1)') me
         else if (me.le.99) then
            write(ppp,'(i2)') me
         else if (me.le.999) then
            write(ppp,'(i3)') me
         else if (me.le.9999) then
            write(ppp,'(i4)') me
         else
            call vberr(' too many processors ')
         end if
         zout = 'out_'//ppp
         if (.not.oroot()) open(6,file=zout)
      end if
c
      write(iwr,601) ipg_nnodes()
      if (zout.eq.'output') then
         write(iwr,602)
      else if (zout.eq.'ibm') then 
         write(iwr,603)
      else
         write(iwr,604)
      end if
         write(iwr,605)
601   format(' *************************************************',
     1      /' **  Parallel Processing on',i6,' processors    **')
602   format(' **  output is done only by the root            **')
603   format(' **  output is done by all (ibm style)          **')
604   format(' **  output is send to file out_process-id      **')
605   format(' *************************************************',/1x)
c
      return
      end
_IFN(ibmmpi)
      subroutine MP_STDOUT_MODE(imode)
      MPIINT imode
      end
_ENDIF
_ENDIF
_IF(atmol)
      subroutine fopen(iunit,name,status,form)
c
      implicit REAL (a-h,o-z), integer (i-n)
c
c.... open a fortran file (may be for paralllel)
c
      common/charsort/sortname(2),scratch,zout
      character*8 zout
      character*44 sortname,scratch
      character*(*) name,status,form
      character*44 filen
      character*4 yyy
c
      filen = name(1:len(name))
_IF(parallel)
      if (status.eq.'priv') then
         n = ipg_nodeid()
         write(yyy,'(i4)') n
         if (n.gt.999) then
            filen = 'N'//yyy//'.'//filen
         else if (n.gt.99) then
            filen = 'N'//yyy(2:4)//'.'//filen
         else if (n.gt.9) then
            filen = 'N'//yyy(3:4)//'.'//filen
         else
            filen = 'N'//yyy(4:4)//'.'//filen
         end if
      end if
_ENDIF
         call strtrm(filen,ll)
         open(iunit,file=filen(1:ll),access='sequential',
     1        form=form)
         rewind iunit
c
      return
      end
_ENDIF
_IF(parallel)
      subroutine pg_dgop_large(type,q,n,op)
c
      parameter (maxbuf = 10000)
c
      integer type, n
      REAL q(n)
      character*(*) op
c
      integer nchunks, nleft, i
c
      nchunks = n / maxbuf
      nleft    = mod(n,maxbuf)
c     
      do i=0,nchunks-1
        call pg_dgop(type,q(i*maxbuf+1),maxbuf,op)
      end do 
      call pg_dgop(type,q(nchunks*maxbuf+1),nleft,op)
      return
      end 
_IF(peigss)
      subroutine cleargas()
c
      logical oga, pg_destroy
INCLUDE(common/gainfo)

      oga = pg_destroy(iga_h)
      if(.not.oga) call ga_error
     &                 ('GA: failed to destroy hmat_ga',iga_h)
      oga = pg_destroy(iga_s)
      if (.not.oga) call ga_error
     &                  ('GA: failed to destroy smat_ga',iga_s)
      return
      end 
_ENDIF      
_ENDIF
c***********************************************************************
      subroutine prweig(weight,irumer,ilead,nstruc,weightgn)
c
      implicit none
c
      REAL weight(*),weightgn(*)
      integer irumer(*),nstruc,ilead(*)
c
      common /indep/ time0,nelec,mult
      REAL time0
      integer nelec,mult
      integer kscra7vb,k7bond
INCLUDE(../m4/common/iofile)
INCLUDE(common/turtleparam)
      REAL prlev
      common/leading/prlev
c
      integer i,j,nbonds
c...   nbonds is strictly local / not related to /vbbonds/
      character*1 ab(2)
      data ab/'a','b'/
c
INCLUDE(common/logics)
c
      if (.not.manual) then
c     do not print special list with manual configurations
       nbonds=(nelec-mult+1)/2
       k7bond = kscra7vb('k7bond',nstruc*nelec,'i','r')
       call readis(irumer,nstruc*nelec,num8)
       if (nstruc.gt.100) prlev = 0.005
       write(iwr,2) prlev
  2    format(/,' list contains diagrams for G-N weights larger than ',
     &        f11.8)
       if (allrum) then
         write(iwr,5) 
  5      format(/,'  s#  |weight Ch-C |weight G-N |rumer diagrams',/,
     &    '------+------------+-----------+--',
     &    '------------------------------------------',
     &    '-------------------------------')
          do i=1,nstruc
            if ( weightgn(i) .gt. prlev) then
              write(iwr,10) i,weight(i),weightgn(i),
     &        ('(',abs(irumer((i-1)*nelec+2*
     &        j-1)),' ',' ',abs(irumer((i-1)*nelec+2*j)),')',j=1,
     &        nbonds),(' ',irumer((i-1)*nelec+j),' ',
     & j=2*nbonds+1,nelec)
            endif
 10         format(i5,'  ',f10.6,'  ',f10.6,(t34,22(a1,i2,a1)))
          enddo
       else
         call leadsup(irumer,ilead,nstruc*nelec)
         write(iwr,15) 
 15      format(/,' s# |weight Ch-C |weight G-N |branching diagrams',/,
     &    '----+------------+-----------',
     &    '+--------------------------------------------',
     &    '-----------------------------------')
         do i=1,nstruc
           if (weightgn(i) .gt. prlev) then
             write(iwr,20) i,weight(i),weightgn(i),
     &       (ab(ilead((i-1)*nelec+j)),j=1,
     &       nelec)
           endif
  20       format (i3,'  ',f10.6,'  ',f10.6,t32,100a1)
         enddo
       endif
      endif
      return
      end
      integer function kscra7vb(name,nw,type,rw)
c
      implicit none
c
c...  dumpfile like organisation of vb file - scra7vb
c
c...  type : 'i' or 'I' : integer
c...         'r' or 'R' : real*8
c...         'b' or 'B' : blocks
c...  rw  : intention
c...        i : initialise system / set kbl7vb to nw / no search
c...        j : set kbl7vb to nw / no search
c...        k : return the current value of kbl7vb no search
c...        w : write
c...        r : read
c...        p : (postpone) - w do not know the nw yet; becomes clear on 'a'
c...        a : finishes 'p' - name must be the same; no kscra7vb calls in between
c...        n : return start of name, without size check ; no search
c...        c : return block beyond last block in nw ; no search
c...            This could be used for more precise checking, but requires precise dims
c
c
INCLUDE(../m4/common/iofile)
INCLUDE(common/scra7vb)
      integer nkkk,flip_scra7vb,iposun,kbl7vb
      parameter (nkkk=30)
      integer nw,nword,k,l,lll,lensec,nipw,npoint
      integer kpoint(nkkk),length(nkkk),size(nkkk)
      character*(*) name,name2
      character*30 prep,point(nkkk)
      character*1 type,rw
      save point,prep,length,size,kpoint,npoint,kbl7vb
c
      kscra7vb = 0
      if (rw.eq.'i') then
         do k=1,nkkk
            point(k) = ' '
            length(k) = 0
            size(k) = 0
            kpoint(k) = 0
         end do
         npoint = 0
         kbl7vb = nw
         return
      else if (rw.eq.'j') then
         kbl7vb = nw
         return
      else if (rw.eq.'k') then
         kscra7vb = kbl7vb
         return
      end if
c
      if (type.eq.'i'.or.type.eq.'I') then
         nword = (nw-1)/nipw() + 1
      else if (type.eq.'b'.or.type.eq.'B') then
         nword = nw * 511
         if (lensec(nword).ne.nw) call caserr('kscravb lensec problem')
      else
         nword = nw
      end if
c
      do k=1,npoint
         if (point(k).eq.name) go to 20
      end do
c
      if (rw.eq.'c') then
c...   not found - return -1 -1
         nw= -1
         kscra7vb = -1
         return
      else if (rw.eq.'n') then
c...   not found - return -1 
         kscra7vb = -1
         return
      end if
c
c...  new point
c
      npoint = npoint + 1
      k = npoint
      if (k.gt.nkkk) then
         write(iwr,*) 'kscra7vb overflow on adding ',name
         call caserr('kscra7vb overflow')
      end if
      point(k) = name
10    kpoint(k) = kbl7vb
      if (rw.eq.'p') then
c...    prepare action ; a finishes ; always make a new section
        prep = name
        go to 30
      else if (rw.eq.'a') then
         if (name.ne.prep) call caserr('prepname wrong')
         length(k) = -1
         kbl7vb = iposun(num8) 
         size(k) = kbl7vb - kpoint(k)
         prep = ' '
         goto 35
      else if (rw.eq.'r') then
         write(iwr,*) 'read before write in kscra7vb on ',name
         call caserr('read before write in ksra7vb')
      else if (rw.eq.'w') then
         length(k) = nword
         size(k) = lensec(nword)
         kbl7vb = kbl7vb + size(k)
         go to 30
      else
         call caserr(' illegal rw in kscra7vb')
      end if
c
c...  existing point
c
20    if (rw.eq.'a') then
         go to 10
      else if  (rw.eq.'p') then
         if (prep.eq.' ') go to 10
         write(iwr,*) ' p (',prep,') should be closed by a before ',name
         call caserr('p not finished  kscra7vb')
      else if (rw.eq.'c') then
         nw = kpoint(k) + size(k)
         go to 35
      else if (rw.eq.'n') then
         go to 35
      else if (rw.eq.'r') then
         if (nword.ne.length(k).and.length(k).ne.-1) then
            write(iwr,*) ' ** read kscravb section length ',nword,
     1                   ' not ',length(k)
            call caserr('trying to read more or less in kscra7vb')
         end if
         lll = lensec(nword)
         if (lll.gt.size(k)) then
            write(iwr,*) ' ** kscravb section from ',size(k),' to ',lll
            call caserr('trying to read more in kscra7vb')
         end if
      else if (rw.eq.'w') then
         lll = lensec(nword)
         if (lll.gt.size(k)) then
            write(iwr,*) ' ** kscravb section from ',size(k),' to ',lll
            go to 10
         end if
         length(k) = nword
      else 
         call caserr('illegal rw in exist scra7vb')
      end if
c
30    call search(kpoint(k),num8)
35    kscra7vb = kpoint(k)
c
      return
c
      entry flip_scra7vb(name,name2)
c...  flip the two
      do k=1,npoint
         if (point(k).eq.name) go to 80
      end do
      call caserr(' name not found in nflip_scra7')
80    do l=1,npoint
         if (point(l).eq.name2) go to 90
      end do
      call caserr(' name2 not found in nflip_scra7')
90    lll = kpoint(l)
      kpoint(l) = kpoint(k)
      kpoint(k) = lll
      if (length(k).ne.length(l)) call caserr('flip_scr7vb error')
c
      return
      end         
