dnl this is the machine dependent file for DEC alpha AXP osf workstation
define(GEN_OPTIONS,`dec,osf,alpha,blas,cio')dnl
define(GEN_MACHINE,`d')dnl
define(REAL,real*8)dnl
define(COMPLEX,complex*16)dnl
changequote(!,!)dnl
define(INCLUDE,!      include '$1'!)dnl
changequote()dnl
define(dfloat,dble($1))dnl
