dnl this is the machine dependent file for ksr
define(GEN_OPTIONS,`ksr,cio,veclib,blas')dnl
define(GEN_MACHINE,`k')dnl
define(REAL,real*8)dnl
define(INTEGER,integer)dnl
define(COMPLEX,complex*16)dnl
changequote(!,!)dnl
define(INCLUDE,!      include '$1'!)dnl
changequote()dnl
define(idamax,isamax($@))dnl
define(ddot,sdot($@))dnl
define(dscal,sscal($@))dnl
define(dswap,sswap($@))dnl
define(dgemm,sgemm($@))dnl
define(dgemv,sgemv($@))dnl
define(dnrm2,snrm2($@))dnl
define(dasum,sasum($@))dnl
define(drot,srot($@))dnl
define(daxpy,saxpy($@))dnl
define(dcopy,scopy($@))dnl
dnl --------------------------------------------------
dnl   functions switches which are conditional on the 
dnl   stride arguments being 1 - the stride args are 
dnl   not 1 in the fortran subroutine declaration 
dnl   so they survives
dnl --------------------------------------------------
define(dand,and($@))dnl
define(dor, or($@) )dnl
define(dxor,xor($@))dnl
define(shiftr,rshift($@))dnl
define(shiftl,lshift($@))dnl
