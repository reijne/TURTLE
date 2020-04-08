dnl this is the machine dependent file for cyber 205
define(GEN_OPTIONS,`cyber205')dnl
define(GEN_MACHINE,`u')dnl
define(REAL,real)dnl
define(COMPLEX,complex)dnl
changequote(!,!)dnl
define(INCLUDE,!      include '$1'!)dnl
changequote()dnl
dnl ----- the following are fortran intrinsics -------
dnl       precision change only
define(dlog,alog($@))dnl
define(dlog10,alog10($@))dnl
define(dexp,exp($@))dnl
define(dsqrt,sqrt($@))dnl
define(datan,atan($@))dnl
define(datan2,atan2($@))dnl
define(dsin,sin($@))dnl
define(dcos,cos($@))dnl
define(dabs,abs($@))dnl
define(dmax1,amax1($@))dnl
define(dmin1,amin1($@))dnl
define(idint,int($@))dnl
define(dmod,amod($@))dnl
define(dfloat,float($@))dnl
define(dsign,sign($@))dnl
define(zexp,cexp($@))dnl
define(zlog,clog($@))dnl
define(zlog10,clog10($@))dnl
define(zsqrt,csqrt($@))dnl
define(zcos,ccos($@))dnl
define(zsin,csin($@))dnl
define(dimag,aimag($@))dnl
define(dtan,tan($@))dnl
define(dasin,asin($@))dnl
define(dacos,acos($@))dnl
dnl ---- blas names ----
define(idmax,ismax($@))dnl
define(idamax,isamax($@))dnl
define(dswap,sswap($@))dnl
define(dsum,ssum($@))dnl
define(dasum,sasum($@))dnl
define(idmin,ismin($@))dnl
define(idamin,isamin($@))dnl
define(ddoti,sdoti($@))dnl
define(drot,srot($@))dnl
dnl ---  blas names removed for cyber ----
dnl define(dnrm2,snrm2($@))dnl
dnl ------- assembler --------
define(vclr,szero($1,$3))dnl
dnl ------ cray/cyber with exact arg match
define(dgthr,gather($1,$3,$2,$4))dnl
define(dsctr,scatter($1,$4,$3,$2))dnl
define(vfill,sfill($4,$1,$2,$3))dnl
dnl --------------------------------------------------
dnl   functions switches which are conditional on the 
dnl   stride arguments being 1 - the stride args are 
dnl   not 1 in the fortran subroutine declaration 
dnl   so they survives
dnl --------------------------------------------------
define(dcopy,`ifelse(TEST_ARGS($3,1,$5,1),fmove($2,$4,$1),`scopy($1,$2,$3,$4,$5)')')dnl
define(ddot,`ifelse(TEST_ARGS($3,1,$5,1),vecsum($2,$4,$1),`sdot($1,$2,$3,$4,$5)')')dnl
define(daxpy,`ifelse(TEST_ARGS($4,1,$6,1),triad($1,$2,$5,$3),`saxpy($1,$2,$3,$4,$5,$6)')')dnl
define(dscal,`ifelse(TEST_ARGS($4,1),scaler($1,$2,$3,$3),`sscal($1,$2,$3,$4)')')dnl
define(vmul,`ifelse(TEST_ARGS($2,1,$4,1,$6,1),vvtv($7,$5,$1,$3),``vmul($1,$2,$3,$4,$5,$6,$7)'')')dnl
define(vsma,`ifelse(TEST_ARGS($2,1,$5,1,$7,1),gtriad($8,$3,$6,$4,$1),``vsma($1,$2,$3,$4,$5,$6,$7,$8)'')')dnl
define(vadd,`ifelse(TEST_ARGS($2,1,$4,1,$6,1),addvec($5,$1,$3,$7),``vadd($1,$2,$3,$4,$5,$6,$7)'')')dnl
define(dnrm2,`ifelse(TEST_ARGS($3,1),sqrt(vecsum($2,$2,$1)),`snrm2($1,$2,$3)')')dnl
