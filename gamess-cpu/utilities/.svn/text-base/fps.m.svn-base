dnl this is the machine dependent file for fps 164
define(GEN_OPTIONS,`fps,164,vaxh')dnl
define(GEN_MACHINE,`f')dnl
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
define(dsqrt,sqrt($1))dnl
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
define(ddot,sdot($@))dnl
define(dscal,sscal($@))dnl
define(dswap,sswap($@))dnl
define(dnrm2,snrm2($@))dnl
define(dsum,ssum($@))dnl
define(dasum,sasum($@))dnl
define(idmin,ismin($@))dnl
define(idamin,isamin($@))dnl
define(ddoti,sdoti($@))dnl
define(drot,srot($@))dnl
define(daxpy,saxpy($@))dnl
define(dcopy,scopy($@))dnl
define(dand,and($@))dnl
define(dor, or($@) )dnl
define(dxor,xor($@))dnl
dnl ------- names specific to the fps ----
define(setsto,vfill($2,$3,1,$1))dnl
define(dgthr,gather($1,$3,$2,$4))dnl
define(dsctr,scatter($1,$4,$3,$2))dnl
define(trnsps,mtrans($1,1,$2,1,$4,$3))dnl
dnl ------ fps quirk
define(vsub,``vsub($3,$2,$1,$4,$5,$6,$7)'')dnl
