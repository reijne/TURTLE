dnl A macro to determine which compiler is being used, in order that
dnl different flags can be set

AC_DEFUN([TW_FC_ID], [
AC_REQUIRE([AC_PROG_FC])

FC_ID=

dnl Firstly go by compiler name.

case $FC in 
   
   g77*)
      FC_ID=G77
      ;;

   gfortran*)
      FC_ID=Gfortran
      ;;

   if*)
      FC_ID=Intel
      ;;

   lf9*)
      FC_ID=Lahey
      ;;
   
   pgf*)
      FC_ID=Portland
      ;;

   xlf*)
      FC_ID=Xlf 

esac

dnl then try and disambiguate all f77, f90, and f95 types.
dnl We should have a choice between
dnl nag. absoft. sun. sgi. digital. hp. cray. ...?

if test x$FC_ID = x; then
   tw_fc_v_output=$($FC -V 2>&1 )
   if test $?; then
      case $tw_fc_v_output in
         *NAG*)
            FC_ID=Nag
            ;;
         *Sun*)
            FC_ID=Sun # there's more than one compiler here ...
            ;;
      esac
   fi
fi
 if test x$FC_ID = x; then
   tw_fc_v_output=$($FC -version 2>&1)
   if test $?; then
      case $tw_fc_v_output in
         *Compaq*)
            FC_ID=Digital
            ;;
         *Digital*)
            FC_ID=Digital
            ;;
         *SGI*)
            FC_ID=SGI
            ;;
      esac
   fi
fi   
   
AS_IF([test x$FC_ID != x],
      [AC_MSG_NOTICE([$FC seems to be a $FC_ID compiler])],
      [FC_ID=unknown; AC_MSG_NOTICE([Could not determine type of compiler])])

dnl for more fun, try and get the version number now ...


])# TW_FC_ID
