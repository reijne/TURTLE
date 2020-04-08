      subroutine missing_module(name)
      implicit none
c
      integer ERR_NO_CODE
      parameter(ERR_NO_CODE=0)
c
c
c ==== if the error will occur on all nodes, we can ====
c      handle the output more cleanly
c
      integer ERR_SYNC, ERR_ASYNC
      parameter(ERR_ASYNC=0)
      parameter(ERR_SYNC=301)
c
c ==== legal error classes - these control an explanatory  ====
c      message (one per class, see machscf.m).
c
      integer ERR_NO_CLASS, ERR_INCOMPREHENSIBLE,
     &  ERR_INCOMPATIBLE, ERR_DIMENSION, ERR_FIXED_DIMENSION,
     &  ERR_UNIMPLEMENTED, ERR_INTERNAL,  ERR_USER_DIMENSION,
     &  ERR_SWITCHED_OUT, ERR_UNLUCKY

      parameter(ERR_NO_CLASS=0)
c - fix input errors and try again
      parameter(ERR_INCOMPREHENSIBLE=101)
c - incompatible 
      parameter(ERR_INCOMPATIBLE=102)
c - code dimensions exceeded - redimension code
      parameter(ERR_DIMENSION=103)
c - code dimensions exceeded - redimension or contact support
      parameter(ERR_FIXED_DIMENSION=104)
c - request unimplemented option
      parameter(ERR_UNIMPLEMENTED=105)
c - internal error - please report to support
      parameter(ERR_INTERNAL=106)

c - input dimension exceeded - change allocation
c   in input file and re-run
      parameter(ERR_USER_DIMENSION=107)

c - option disabled at configure stage
      parameter(ERR_SWITCHED_OUT=108)

c - some algorithmic or case-specififc problem
      parameter(ERR_UNLUCKY=109)

c
c  System message flag
c
      integer ERR_NO_SYS, ERR_SYS
      parameter(ERR_NO_SYS=0)
      parameter(ERR_SYS=201)
      character errstr*132, name*(*)

      integer lstchr
      external lstchr

      errstr = "Feature "//name(1:lstchr(name))//
     & " is not available in this build"

      call gamerr(errstr,
     &     ERR_NO_CODE, ERR_SWITCHED_OUT, ERR_SYNC, ERR_NO_SYS)

      end

c
c  The first set are excluded from all builds
c
c Plot option (convex specific)

      subroutine gamplt
      entry  iniplt
      entry stoplt
      entry cntour
      entry arch3d
      call missing_module('plot module')
      return
      end

c MOPAC
      subroutine mopin
      entry mopac
      call missing_module('MOPAC')
      return
      end




      subroutine masscf(core,total_energy)
      entry masscf_input
      call missing_module('MASSCF')
      end

      subroutine ver_masscf(s,r,d)
      character s*(*),r*(*),d*(*)
      s=" "
      end







