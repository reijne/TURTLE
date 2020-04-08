      subroutine missing_module(name)
      implicit none
INCLUDE(common/errcodes)
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
_IFN(mopac)
      entry mopac
      call missing_module('MOPAC')
_ENDIF
      return
      end

_IFN(drf)
      subroutine mcx
      call missing_module('DRF')
      end
_ENDIF

_IF(qmmm)
      subroutine ircdat
      entry irc
      call missing_module('IRC')
      end
      subroutine optjs
      call missing_module('JS Optimisation')
      end
      subroutine ver_optef(s,r,d)
      character s*(*), r*(*), d*(*)
      s=" "
      end
_ENDIF

_IFN(vb)
      subroutine vbstart(core,n)
      entry tranvbini        
      call missing_module('VB')
      end
      subroutine servec
      call missing_module('VB')
      end
_ENDIF

_IFN(masscf)
      subroutine masscf(core,total_energy)
      entry masscf_input
      call missing_module('MASSCF')
      end

      subroutine ver_masscf(s,r,d)
      character s*(*),r*(*),d*(*)
      s=" "
      end
_ENDIF

_IFN(zora)
      subroutine zora
      entry scf_so
      entry basis_zora
      entry scale_z
      entry scale_uhf
      entry coul_z
      entry swap_zora
      entry zora_int
      call missing_module('ZORA' )
      end

      subroutine ver_zora(s,r,d)
      character s*(*),r*(*),d*(*)
      s=" "
      end
_ENDIF

_IF(ga_mp2_build,ga_ci_build)
      subroutine ver_mp2a(s,r,d)
      character s*(*), r*(*), d*(*)
      s=" "
      end
      subroutine ver_mp2b(s,r,d)
      character s*(*), r*(*), d*(*)
      s=" "
      end
      subroutine ver_mp2c(s,r,d)
      character s*(*), r*(*), d*(*)
      s=" "
      end
      subroutine ver_mp2d(s,r,d)
      character s*(*), r*(*), d*(*)
      s=" "
      end
      subroutine ver_mp3(s,r,d)
      character s*(*), r*(*), d*(*)
      s=" "
      end
      subroutine ver_secmp2(s,r,d)
      character s*(*), r*(*), d*(*)
      s=" "
      end

      subroutine mp2eng
      entry mp2ddrv
      entry mp2dmd
      entry mpdmdw
      call missing_module('MP2 second derivatives')
      end

_IF(qmmm)
c
c specific exclusion of rpa and new direct scf
c codes in qm/mm case 
c
      subroutine rpinit
      entry rpdriv
      call missing_module('RPA')
      end
      subroutine ver_rpa(s,r,d)
      character s*(*), r*(*), d*(*)
      s=" "
      end
      subroutine ver_dirrpa(s,r,d)
      character s*(*), r*(*), d*(*)
      s=" "
      end
      subroutine drcdrv
      call missing_module('New Direct')
      end
      subroutine ver_direct(s,r,d)
      character s*(*), r*(*), d*(*)
      s=" "
      end
_ENDIF
_ENDIF


_IF(base_build,mp2_build,ga_mp2_build,ci_build,ga_ci_build)
c
c  MP2 mode (serial platforms)
c
c
c CCSD
c
      subroutine ccsdin
      entry titandrv
      call missing_module('CCSD' )
      end

      subroutine ver_ccsd(s,r,d)
      character s*(*), r*(*), d*(*)
      s=" "
      end
      subroutine ver_tsortc(s,ls,r,lr,d,ld)
      character s*(*), r*(*), d*(*)
      s=" "
      end
      subroutine ver_tsort(s,r,d)
      character s*(*), r*(*), d*(*)
      s=" "
      end
      subroutine ver_nvccsd(s,r,d)
      character s*(*), r*(*), d*(*)
      s=" "
      end

c
c Model (AMBER QM/MM)
c
      subroutine model
      call missing_module('Model QM/MM' )
      end
      subroutine ver_model(s,r,d)
      character s*(*), r*(*), d*(*)
      s=" "
      end
c
c MCLR/RPA
c
_IF(ga_mp2_build,ga_ci_build)
      subroutine lrinit
      entry lrmain
      call missing_module('MCLR')
      end

      subroutine ver_mclr(s,r,d)
      character s*(*), r*(*), d*(*)
      s=" "
      end
_ELSE
      subroutine lrinit
      entry rpinit
      entry rpdriv
      entry lrmain
      call missing_module('MCLR')
      end

      subroutine ver_rpa(s,r,d)
      character s*(*), r*(*), d*(*)
      s=" "
      end
      subroutine ver_dirrpa(s,r,d)
      character s*(*), r*(*), d*(*)
      s=" "
      end
      subroutine ver_mclr(s,r,d)
      character s*(*), r*(*), d*(*)
      s=" "
      end
_ENDIF

c
c Greens function methods
c
_IF(ga_ci_build)
      subroutine tda
      call missing_module('Greens Function')
      end

      subroutine ver_tdaf(s,r,d)
      character s*(*), r*(*), d*(*)
      s=" "
      end
_ELSE
      subroutine inip
      entry tda
      entry gf
      call missing_module('Greens Function')
      end

      subroutine ver_tdaf(s,r,d)
      character s*(*), r*(*), d*(*)
      s=" "
      end
      subroutine ver_gff(s,r,d)
      character s*(*), r*(*), d*(*)
      s=" "
      end
_ENDIF
c
c Full-CI
c
      subroutine fci
      entry fulcin
      call missing_module('Full Ci')
      end

      subroutine ver_fullci(s,r,d)
      character s*(*), r*(*), d*(*)
      s=" "
      end
_IFN(ga_mp2_build,ga_ci_build)
c
c New Direct-SCF
c
      subroutine drcdrv
      call missing_module('New Direct')
      end

      subroutine ver_direct(s,r,d)
      character s*(*), r*(*), d*(*)
      s=" "
      end
_ENDIF

_IFN(vb)
c
c CASSCF module
c
      subroutine casscf
      call missing_module('CASSCF')
      end
      subroutine ver_casb(s,r,d)
      character s*(*), r*(*), d*(*)
      s=" "
      end
      subroutine ver_casa(s,r,d)
      character s*(*), r*(*), d*(*)
      s=" "
      end

_ENDIF

c
c MSCSF
c
      subroutine mcstar
      entry mcdump
      entry mcgrad
      entry mcanal
      entry wvfn
_IFN(charmm)
c sort clashes with a charmm utility
      entry sort
_ENDIF
      entry redun
      entry fmtp
      entry dmpini
      entry mcprorb
      call missing_module('MCSCF')
      end

      subroutine ver_mcscfa(s,r,d)
      character s*(*), r*(*), d*(*)
      s=" "
      end
      subroutine ver_mcscfc(s,r,d)
      character s*(*), r*(*), d*(*)
      s=" "
      end
      subroutine ver_mcscfb(s,r,d)
      character s*(*), r*(*), d*(*)
      s=" "
      end

c
c Saunders/van Lenthe Direct-CI
c
_IFN(ga_ci_build)
      subroutine cidata
      entry cizero
_IF(charmm)
      entry omitdirect      
_ELSE
      entry direct
_ENDIF
      call missing_module('Direct-CI')
      end
      subroutine ver_dircta(s,r,d)
      character s*(*), r*(*), d*(*)
      s=" "
      end
      subroutine ver_dirctb(s,r,d)
      character s*(*), r*(*), d*(*)
      s=" "
      end
      subroutine ver_dirctc(s,r,d)
      character s*(*), r*(*), d*(*)
      s=" "
      end
      subroutine ver_dirctd(s,r,d)
      character s*(*), r*(*), d*(*)
      s=" "
      end
_ENDIF

c
c MRDCI modules
c 

      subroutine pmrdci
      entry nmrdci
      entry dmrdci
      entry smrdci
      entry tdata
      entry tmrdci
      entry amrdci
      entry moment
      entry tabci
      call missing_module('MRDCI')
      end

      subroutine ver_mrdci7(s,r,d)
      character s*(*), r*(*), d*(*)
      s=" "
      end
      subroutine ver_mrdci6(s,r,d)
      character s*(*), r*(*), d*(*)
      s=" "
      end
      subroutine ver_mrdci5(s,r,d)
      character s*(*), r*(*), d*(*)
      s=" "
      end
      subroutine ver_mrdci4(s,r,d)
      character s*(*), r*(*), d*(*)
      s=" "
      end
      subroutine ver_mrdci3(s,r,d)
      character s*(*), r*(*), d*(*)
      s=" "
      end
      subroutine ver_mrdci2(s,r,d)
      character s*(*), r*(*), d*(*)
      s=" "
      end
      subroutine ver_mrdci1(s,r,d)
      character s*(*), r*(*), d*(*)
      s=" "
      end

c New MRDCI

_IFN(mrdci)
      subroutine mrdcin2
      entry mrdcim2
      entry parkinc
      entry sortdriver
      call missing_module('new MRDCI')
      end
_ENDIF

_ENDIF

_IF(base_build)
c
c Stubs for code that is present in the MP2 build
c Including inpci
c

      subroutine multi
      entry grmp23
      entry emp23
      entry mp2eng
      entry uhfmp2
      entry uhfmp3
      entry mptran
      entry adapti
_IFN(vb)
      entry inporb
      entry closbf
      entry closbf2
_ENDIF
      entry raman
      entry hfdder
      entry hfhpol
      entry hfpder
      entry chfidr
      entry mp2dd
      entry anairr
      entry scfdd
      entry chfndr
      entry trnfkd
      entry indxsy
      entry mpdder
      entry trndrv
      entry incas1
      entry ciin
      entry multin
      entry mrdcin
_IFN(vb)
      entry setbfb
      entry setbfa
_ENDIF
      entry incas3
      entry incas2
      entry mrdcim
      call missing_module('MP2/SecDer')
      end

_IFN(vb)
      subroutine ver_util4(s,r,d)
      character s*(*), r*(*), d*(*)
      s = " "
      end

      subroutine ver_util5(s,r,d)
      character s*(*), r*(*), d*(*)
      s = " "
      end
_ENDIF

      subroutine ver_util6(s,r,d)
      character s*(*), r*(*), d*(*)
      s = " "
      end

      subroutine ver_util7(s,r,d)
      character s*(*), r*(*), d*(*)
      s = " "
      end

      subroutine ver_util8(s,r,d)
      character s*(*), r*(*), d*(*)
      s = " "
      end


      subroutine ver_tran4(s,r,d)
      character s*(*), r*(*), d*(*)
      s = " "
      end


      subroutine ver_sec1e(s,r,d)
      character s*(*), r*(*), d*(*)
      s = " "
      end

      subroutine ver_sec2e(s,r,d)
      character s*(*), r*(*), d*(*)
      s = " "
      end

      subroutine ver_secmp2(s,r,d)
      character s*(*), r*(*), d*(*)
      s = " "
      end


      subroutine ver_mp3(s,r,d)
      character s*(*), r*(*), d*(*)
      s = " "
      end

      subroutine ver_mp2d(s,r,d)
      character s*(*), r*(*), d*(*)
      s = " "
      end

      subroutine ver_mp2c(s,r,d)
      character s*(*), r*(*), d*(*)
      s = " "
      end

      subroutine ver_mp2b(s,r,d)
      character s*(*), r*(*), d*(*)
      s = " "
      end

      subroutine ver_mp2a(s,r,d)
      character s*(*), r*(*), d*(*)
      s = " "
      end


_IFN(vb)
      subroutine ver_mainci(s,r,d)
      character s*(*), r*(*), d*(*)
      s = " "
      end

      subroutine ver_machci(s,r,d)
      character s*(*), r*(*), d*(*)
      s = " "
      end
_ENDIF


      subroutine ver_intege(s,r,d)
      character s*(*), r*(*), d*(*)
      s = " "
      end

_IFN(drf,vb)
cjmht drf can be built with base build but it requires qmint 
c     and qmsint from integd so this file is included in the build
      subroutine ver_integd(s,r,d)
      character s*(*), r*(*), d*(*)
      s = " "
      end
_ENDIF
      subroutine ver_inpci(s,r,d)
      character s*(*), r*(*), d*(*)
      s = " "
      end

      subroutine ver_index4(s,r,d)
      character s*(*), r*(*), d*(*)
      s = " "
      end

      subroutine ver_drvmp(s,r,d)
      character s*(*), r*(*), d*(*)
      s = " "
      end

      subroutine ver_derdrv(s,r,d)
      character s*(*), r*(*), d*(*)
      s = " "
      end

      subroutine ver_cphf(s,r,d)
      character s*(*), r*(*), d*(*)
      s = " "
      end

      subroutine ver_anald(s,r,d)
      character s*(*), r*(*), d*(*)
      s = " "
      end

      subroutine ver_anale(s,r,d)
      character s*(*), r*(*), d*(*)
      s = " "
      end

_ENDIF


_IF(t3d)
c
c Machine specific dummy routines:
c
      subroutine random
      call caserr('random')
      end
_ENDIF
