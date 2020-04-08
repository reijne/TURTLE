dnl this is the machine dependent file for iris with mopac
define(GEN_OPTIONS,`sgi,cio,mopac')dnl
define(GEN_MACHINE,`g')dnl
define(REAL,real*8)dnl
define(COMPLEX,complex*16)dnl
changequote(!,!)dnl
define(INCLUDE,!      include '$1'!)dnl
changequote()dnl
define(dfloat,dble($1))dnl
