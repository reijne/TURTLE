/*
 $Id: globalp.c.h,v 1.2 2000-10-26 15:38:13 psh Exp $

 file globalp.c.h 
 */

/* GAMESS-UK - following if test added */

#ifndef MA_DEFINES_TYPES

#ifdef  STD_INT
typedef int    Integer;
typedef unsigned int unInteger;
#else
typedef long   Integer;
typedef unsigned long unInteger;
#endif

#ifdef  STD_DBL
typedef   double         DoublePrecision;
#else
typedef   long double    DoublePrecision;
#endif

#endif

#include "blas_lapack.h"
#include "peigs_types.h"

