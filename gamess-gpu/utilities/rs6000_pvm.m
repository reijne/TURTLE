dnl this is the machine dependent file for the rs6000, with PVM communications
define(GEN_OPTIONS,`rs6000,cio,parallel,pvm')dnl
define(GEN_MACHINE,`r')dnl
define(REAL,real*8)dnl
define(COMPLEX,complex*16)dnl
changequote(!,!)dnl
define(INCLUDE,!      include '$1'!)dnl
changequote()dnl
define(dfloat,dble($1))dnl
define(SEND,sndpvm(`$1',`$2',`$3',`$4',`$5'))dnl
define(RECV,rcvpvm(`$1',`$2',`$3',`$4',`$5',`$6',`$7'))dnl
