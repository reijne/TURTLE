# _AC_FC_CHECK_SIZEOF (TYPE, [IGNORED], [INCLUDES = DEFAULT-INCLUDES])
# --------------------------------------------------------------------
#
# Determine the size of types and the return of the len intrinsic for Fortran.
# First the binary representation of TYPE is written to conftest.unf between
# two a. After that conftest.unf is reread and the distance between 2 'a' gives
# the type length. In the end the detected size is put in conftest.res.
#
# It will fail if there is a 'a' in binary format in the header of 
# conftest.unf, if there is a 'a' or a newline in the binary representation.

AC_DEFUN([_AC_FC_CHECK_SIZEOF],
[_AC_FORTRAN_ASSERT()dnl
AC_CACHE_CHECK([for $[]_AC_FC[] size of $1],
               [ac_cv_[]_AC_LANG_ABBREV[]_sizeof_$1],
[AC_RUN_IFELSE([AC_LANG_SOURCE([[      program typesize]

      $1 [ttype
      integer j,k,l,stat,length(2)
      character string*30
      logical intype
      ttype = 1
      l = 1

      open (10,FORM='UNFORMATTED',FILE='conftest.unf',STATUS='NEW',
     +IOSTAT=stat,ERR=100)
      if (stat .ne. 0) goto 100
      write(10) 'a',ttype,'a','a',len(string),'a'
      close(10,IOSTAT=stat,ERR=100)
      if (stat .ne. 0) goto 100

1000  format (A)
      open (10,FILE='conftest.unf',STATUS='OLD',IOSTAT=stat,ERR=100)
      if (stat .ne. 0) goto 100
      length = 0
      intype = .false.
      do 20 k=1,6
      read (10, 1000, END=110) string
110   continue
      do 10 j=1,len(string)
      if (intype) then
      if (string(j:j) .eq. 'a') then
      if (l .eq. 2) then
      goto 120
      else
      l = l + 1
      intype = .false.
      endif
      else
      length(l) = length(l) + 1
      endif
      else
      if (string(j:j) .eq. 'a') then
      intype = .true.
      endif
      endif 
10    continue
20    continue
120   close(10,IOSTAT=stat,ERR=100)
      if (stat .ne. 0) goto 100

1100  format (I4)

      if (length(1) .gt. 0) then
      open (10,FILE='conftest.res',STATUS='NEW',IOSTAT=stat,ERR=100)
      if (stat .ne. 0) goto 100
      write (10, 1100) length(1)
      close(10)
      endif
      if (length(2) .gt. 0) then
      open (10,FILE='conftest.res2',STATUS='NEW',IOSTAT=stat,ERR=100)
      if (stat .ne. 0) goto 100
      write (10, 1100) length(2)
      close(10)
      endif

100   continue

      stop
      end]
])],
[ if test -f 'conftest.res'; then
   AS_TR_SH(ac_cv_[]_AC_LANG_ABBREV[]_sizeof_$1)=`cat conftest.res`
   rm -f conftest.res
  else AS_TR_SH(ac_cv_[]_AC_LANG_ABBREV[]_sizeof_$1)=0
  fi
  if test -f 'conftest.res2'; then
   AS_TR_SH(ac_cv_[]_AC_LANG_ABBREV[]_sizeof_len)=`cat conftest.res2`
   rm -f conftest.res2
  else AS_TR_SH(ac_cv_[]_AC_LANG_ABBREV[]_sizeof_len)=0
  fi
],
[
  AS_TR_SH(ac_cv_[]_AC_LANG_ABBREV[]_sizeof_$1)=0
  AS_TR_SH(ac_cv_[]_AC_LANG_ABBREV[]_sizeof_len)=0
], [])
rm -f conftest.unf])
   AC_DEFINE_UNQUOTED(AS_TR_CPP([]_AC_FC[]_sizeof_$1), 
     $AS_TR_SH(ac_cv_[]_AC_LANG_ABBREV[]_sizeof_$1),
     [The size of `$1', detected.])
   AC_DEFINE_UNQUOTED(AS_TR_CPP([]_AC_FC[]_sizeof_len), 
     $AS_TR_SH(ac_cv_[]_AC_LANG_ABBREV[]_sizeof_len),
     [The size of `len', detected.])
   AS_TR_SH([]_AC_FC[]_sizeof_$1)=$AS_TR_SH(ac_cv_[]_AC_LANG_ABBREV[]_sizeof_$1)
   AS_TR_SH([]_AC_FC[]_sizeof_len)=$AS_TR_SH(ac_cv_[]_AC_LANG_ABBREV[]_sizeof_len)
])# _AC_FC_CHECK_SIZEOF 


# AC_F77_CHECK_SIZEOF
# -------------------
AC_DEFUN([AC_F77_CHECK_SIZEOF],
[AC_REQUIRE([AC_PROG_F77])dnl
AC_LANG_PUSH(Fortran 77)dnl
_AC_FC_CHECK_SIZEOF([$1],[$2],[$3])
AC_LANG_POP(Fortran 77)dnl
])# AC_F77_CHECK_SIZEOF


# AC_FC_CHECK_SIZEOF
# ------------------
AC_DEFUN([AC_FC_CHECK_SIZEOF],
[AC_REQUIRE([AC_PROG_FC])dnl
AC_LANG_PUSH(Fortran)dnl
_AC_FC_CHECK_SIZEOF([$1],[$2],[$3])
AC_LANG_POP(Fortran)dnl
])# AC_FC_CHECK_SIZEOF
