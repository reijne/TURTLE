dnl this is the machine dependent file for apollo dn10000
define(GEN_OPTIONS,`apollo,cio')dnl
define(GEN_MACHINE,`p')dnl
define(REAL,real*8)dnl
define(COMPLEX,complex*16)dnl
changequote(!,!)dnl
define(INCLUDE,!      include '$1'!)dnl
changequote()dnl
define(dfloat,dble($@))dnl
