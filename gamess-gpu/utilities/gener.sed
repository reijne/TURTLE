#!/bin/sh
#
# script to make generic source from 
# fortran. does not try and line-wrap
#
sed '
:top
   /^c[a-z]*?/{
   b skip2
   }
   /^c/{
   n
   b top
   }
   : skip2
   /^...... *format/{
   n
:cont
   /^..... /{
   b top
   }
   n
   b cont
   }
   s/^      include .sizes./INCLUDE/
   s/^      include .sizes.for./INCLUDE/
   s/^%include .sizes.for./INCLUDE/
   s/real *\* *8 * /REAL  /
   s/real * /REAL  /
   s/it real(/it REAL (/
   s/complex * /COMPLEX      /
   s/\([-+]\{0,1\}[0-9]\{1,\}\.[0-9]\{1,\}\) \{0,\}e \{0,\}\([-+]\{0,1\}[0-9]\{1,\}\)/\1d\2/g
   s/\([-+]\{0,1\}[0-9]\{0,\}\.[0-9]\{1,\}\) \{0,\}e \{0,\}\([-+]\{0,1\}[0-9]\{1,\}\)/\1d\2/g
   s/\([-+]\{0,1\}[0-9]\{1,\}\.[0-9]\{0,\}\) \{0,\}e \{0,\}\([-+]\{0,1\}[0-9]\{1,\}\)/\1d\2/g
s/\([-+ =*/(.,<&]\)alog *(/\1dlog(/g
s/\([-+ =*/(.,<&]\)alog10 *(/\1dlog10(/g
s/\([-+ =*/(.,<&]\)exp *(/\1dexp(/g
s/\([-+ =*/(.,<&]\)sqrt *(/\1dsqrt(/g
s/\([-+ =*/(.,<&]\)atan *(/\1datan(/g
s/\([-+ =*/(.,<&]\)atan2 *(/\1datan2(/g
s/\([-+ =*/(.,<&]\)sin *(/\1dsin(/g
s/\([-+ =*/(.,<&]\)cos *(/\1dcos(/g
s/\([-+ =*/(.,<&]\)abs *(/\1dabs(/g
s/\([-+ =*/(.,<&]\)amax1 *(/\1dmax1(/g
s/\([-+ =*/(.,<&]\)amin1 *(/\1dmin1(/g
s/\([-+ =*/(.,<&]\)int *(/\1idint(/g
s/\([-+ =*/(.,<&]\)amod *(/\1dmod(/g
s/\([-+ =*/(.,<&]\)float *(/\1dfloat(/g
s/\([-+ =*/(.,<&]\)dble *(/\1dfloat(/g
s/\([-+ =*/(.,<&]\)arsin *(/\1dasin(/g
s/\([-+ =*/(.,<&]\)arcos *(/\1dacos(/g
s/\([-+ =*/(.,<&]\)sign *(/\1dsign(/g
s/\([-+ =*/(.,<&]\)ifix *(/\1idint(/g
s/\([-+ =*/(.,<&]\)cexp *(/\1zexp(/g
s/\([-+ =*/(.,<&]\)clog *(/\1zlog(/g
s/\([-+ =*/(.,<&]\)clog10 *(/\1zlog10(/g
s/\([-+ =*/(.,<&]\)csqrt *(/\1zsqrt(/g
s/\([-+ =*/(.,<&]\)ccos *(/\1zcos(/g
s/\([-+ =*/(.,<&]\)csin *(/\1zsin(/g
s/\([-+ =*/(.,<&]\)aimag *(/\1dimag(/g
s/\([-+ =*/(.,<&]\)tan *(/\1dtan(/g
s/\([-+ =*/(.,<&]\)asin *(/\1dasin(/g
s/\([-+ =*/(.,<&]\)acos *(/\1dacos(/g
s/\([-+ =*/(.,<&]\)ismax *(/\1idmax(/g
s/\([-+ =*/(.,<&]\)isamax *(/\1idamax(/g
s/\([-+ =*/(.,<&]\)sdot *(/\1ddot(/g
s/\([-+ =*/(.,<&]\)sscal *(/\1dscal(/g
s/\([-+ =*/(.,<&]\)sswap *(/\1dswap(/g
s/\([-+ =*/(.,<&]\)sgthr *(/\1dgthr(/g
s/\([-+ =*/(.,<&]\)ssctr *(/\1dsctr(/g
s/\([-+ =*/(.,<&]\)snrm2 *(/\1dnrm2(/g
s/\([-+ =*/(.,<&]\)ssum *(/\1dsum(/g
s/\([-+ =*/(.,<&]\)sasum *(/\1dasum(/g
s/\([-+ =*/(.,<&]\)ismin *(/\1idmin(/g
s/\([-+ =*/(.,<&]\)isamin *(/\1idamin(/g
s/\([-+ =*/(.,<&]\)sdoti *(/\1ddoti(/g
s/\([-+ =*/(.,<&]\)srot *(/\1drot(/g
s/\([-+ =*/(.,<&]\)saxpy *(/\1daxpy(/g
s/\([-+ =*/(.,<&]\)scopy *(/\1dcopy(/g
s/\([-+ =*/(.,<&]\)sdoti *(/\1ddoti(/g
