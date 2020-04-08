dnl this is the machine dependent file for unicos
define(GEN_OPTIONS,`cray,xmp,unicos,vector,parallel,tools')dnl
define(GEN_MACHINE,`c')dnl
define(REAL,real)dnl
define(COMPLEX,complex)dnl
changequote(!,!)dnl
define(INCLUDE,!      include '$1'!)dnl
changequote()dnl
define(cdabs,cabs($@))dnl
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
define(dnint,anint($@))dnl
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
define(idmax,ismax($@))dnl
define(idamax,isamax($@))dnl
define(ddot,sdot($@))dnl
define(dscal,sscal($@))dnl
define(dswap,sswap($@))dnl
define(dgthr,gather($1,$3,$2,$4))dnl
define(dsctr,scatter($1,$4,$3,$2))dnl
define(dnrm2,snrm2($@))dnl
define(dsum,ssum($@))dnl
define(dasum,sasum($@))dnl
define(iand,and($@))dnl
define(idmin,ismin($@))dnl
define(idamin,isamin($@))dnl
define(ddoti,sdoti($@))dnl
define(drot,srot($@))dnl
define(daxpy,saxpy($@))dnl
define(dcopy,scopy($@))dnl
define(vclr,szero($1,$3))dnl
define(vfill,sfill($4,$1,$2,$3))dnl
define(dand,and($@))dnl
define(dor, or($@) )dnl
define(dxor,xor($@))dnl
dnl --------------------------------------------------
dnl   functions switches which are conditional on the 
dnl   stride arguments being 1 - the stride args are 
dnl   not 1 in the fortran subroutine declaration 
dnl   so they survives
dnl --------------------------------------------------
define(vmul,`ifelse(TEST_ARGS($2,1,$4,1,$6,1),vvtv($7,$5,$1,$3),``vmul($1,$2,$3,$4,$5,$6,$7)'')')dnl
define(vadd,`ifelse(TEST_ARGS($2,1,$4,1,$6,1),addvec($5,$1,$3,$7),``vadd($1,$2,$3,$4,$5,$6,$7)'')')dnl
define(vsub,`ifelse(TEST_ARGS($2,1,$4,1,$6,1),subvec($5,$1,$3,$7),``vsub($1,$2,$3,$4,$5,$6,$7)'')')dnl
define(vsma,`ifelse(TEST_ARGS($2,1,$5,1,$7,1),gtriad($8,$3,$6,$4,$1),``vsma($1,$2,$3,$4,$5,$6,$7,$8)'')')dnl
