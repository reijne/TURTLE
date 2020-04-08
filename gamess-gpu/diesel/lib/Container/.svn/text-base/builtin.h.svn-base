// This may look like C code, but it is really -*- C++ -*-

/* 
Copyright (C) 1988, 1992 Free Software Foundation
    written by Doug Lea (dl@rocky.oswego.edu)

This file is part of the GNU C++ Library.  This library is free
software; you can redistribute it and/or modify it under the terms of
the GNU Library General Public License as published by the Free
Software Foundation; either version 2 of the License, or (at your
option) any later version.  This library is distributed in the hope
that it will be useful, but WITHOUT ANY WARRANTY; without even the
implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
PURPOSE.  See the GNU Library General Public License for more details.
You should have received a copy of the GNU Library General Public
License along with this library; if not, write to the Free Software
Foundation, 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.
*/

/*
  arithmetic, etc. functions on built in types
*/


#ifndef _builtin_h
#define _builtin_h 1

#include <stddef.h>
#include <math.h>
#include "../../config.h"

#ifndef __GNUC__
#define __attribute__(x)
#endif

typedef void (*one_arg_error_handler_t)(const char*);
typedef void (*two_arg_error_handler_t)(const char*, const char*);

LONG_INT         gcd(LONG_INT, LONG_INT);
LONG_INT         lg(unsigned LONG_INT); 
double       pow(double, LONG_INT);
LONG_INT         pow(LONG_INT, LONG_INT);

extern "C" double       start_timer();
extern "C" double       return_elapsed_time(double last_time = 0.0);

char*        dtoa(double x, char cvt = 'g', INT width = 0, INT prec = 6);

unsigned INT hashpjw(const char*);
//_G_uint32_t multiplicativehash(_G_int32_t);
unsigned INT foldhash(double);

extern void default_one_arg_error_handler(const char*) __attribute__ ((noreturn));
extern void default_two_arg_error_handler(const char*, const char*) __attribute__ ((noreturn));

extern two_arg_error_handler_t lib_error_handler;

extern two_arg_error_handler_t 
       set_lib_error_handler(two_arg_error_handler_t f);


#if !defined(IV)

inline short abs(short arg) 
{
  return (arg < 0)? -arg : arg;
}

inline INT sign(LONG_INT arg)
{
  return (arg == 0) ? 0 : ( (arg > 0) ? 1 : -1 );
}

inline INT sign(double arg)
{
  return (arg == 0.0) ? 0 : ( (arg > 0.0) ? 1 : -1 );
}

inline LONG_INT sqr(LONG_INT arg)
{
  return arg * arg;
}

#if ! _G_MATH_H_INLINES /* hpux and SCO define this in math.h */
inline double sqr(double arg)
{
  return arg * arg;
}
#endif

inline INT even(LONG_INT arg)
{
  return !(arg & 1);
}

inline INT odd(LONG_INT arg)
{
  return (arg & 1);
}

inline LONG_INT lcm(LONG_INT x, LONG_INT y)
{
  return x / gcd(x, y) * y;
}

inline void (setbit)(LONG_INT& x, LONG_INT b)
{
  x |= (1 << b);
}

inline void clearbit(LONG_INT& x, LONG_INT b)
{
  x &= ~(1 << b);
}

inline INT testbit(LONG_INT x, LONG_INT b)
{
  return ((x & (1 << b)) != 0);
}

#endif
#endif
