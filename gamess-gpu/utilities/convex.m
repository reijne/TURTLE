dnl this is the machine dependent file for iris
define(GEN_OPTIONS,`convex,fortio,vector,splitopts')dnl
define(GEN_MACHINE,`x')dnl
define(REAL,`real*8')dnl
define(COMPLEX,complex*16)dnl
changequote(!,!)dnl
define(INCLUDE,!      include '$1'!)dnl
changequote()dnl
define(dand,iand($1,$2))dnl
define(dor,ior($1,$2))dnl
define(dxor,ieor($1,$2))dnl
