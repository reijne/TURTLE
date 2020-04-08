dnl A macro to set various compiler-dependent things that can't be sensibly
dnl deduced.

AC_DEFUN([TW_FC_ID_FLAGS], [
AC_REQUIRE([TW_FC_ID])

case $FC_ID in

  Absoft)
     FFLAGS_DEBUG="-et -g -Rb -Rc -Rp -Rs"
     ;;

  Digital)
     FFLAGS_FAST=-O2
     FFLAGS_DEBUG="-g -Rabc -ei"
     ;;

  G77)
     ;;

  Gfortran)
     ;;
 
  Intel)
     FFLAGS_DEBUG="-C -g -inline_debug_info"
     ;;

  Lahey)
     FFLAGS_DEBUG="--chk aesux --chkglobal -g --trace"
     FFLAGS_FAST="-O --warn --quiet --tpp --ntrace"
     ;;

  Nag)
     # This is a hack - we should test for these next two
     FCFLAGS="$FCFLAGS -mismatch -kind=byte"
     FFLAGS_MPI="-kind=byte -mismatch"
     FFLAGS_DEBUG="-C=all -g -gline -nan"
     DEFS="$DEFS __NAG__"
     SYS=nag
     ;;
  
  Portland)
     FFLAGS_DEBUG="-g -Mbounds"
     FFLAGS_FAST="-fast"
     ;;

  SGI)
     FFLAGS_DEBUG="-g -O0"
     FFLAGS_FAST="-O3 -OPT:Olimit=0"
     ;;

  Sun)
     FFLAGS_DEBUG="-C -g"
     FFLAGS_FAST="-fast"
     ;;

  Xlf)
     FFLAGS_DEBUG="-g -C -qinitauto -qsave -qmaxmem=16000 -qnolm"
     FFLAGS_FAST="-O3 -qarch=auto -qtune=auto -qcache=auto -qnolm"
     SYS=xlf
     ;;

esac

AC_SUBST(FFLAGS_MPI)

])
