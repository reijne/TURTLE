#ifndef _ma_h
#define _ma_h

/* 
 * $Id: ma.h,v 1.1.1.2 2004-05-07 15:39:25 psh Exp $
 */

/* 
 * Private header file containing symbolic constants and type declarations
 * for internal C routines.
 *
 * This file should only be included by internal C files.
 */

#include "macdecls.h"

/**
 ** function types
 **/

extern Boolean MAi_inform_base();
extern void MAi_summarize_allocated_blocks();
#ifdef __crayx1
#undef  _CRAY
#undef CRAY
#endif


#endif /* _ma_h */
