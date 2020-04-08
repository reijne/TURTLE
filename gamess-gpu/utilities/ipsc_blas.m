dnl this is the machine dependent file for ipsc860
define(GEN_OPTIONS,`ipsc,parallel,cio,blas')dnl
define(GEN_MACHINE,`h')dnl
define(REAL,real*8)dnl
define(COMPLEX,complex*16)dnl
changequote(!,!)dnl
define(INCLUDE,!      include '$1'!)dnl
changequote()dnl
define(dgop,ggop(`$1',`$2',`$3',`$4',`$5'))dnl
define(synch,gsync)dnl
