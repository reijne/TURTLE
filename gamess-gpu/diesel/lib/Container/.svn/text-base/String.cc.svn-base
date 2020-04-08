/* 
Copyright (C) 1988 Free Software Foundation
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
  String class implementation
 */

#include "String.h"
#include <ctype.h>
#include <limits.h>
#include <new>
#include "builtin.h"

using namespace std;

#undef OK

void String::error(const char* msg) const
{
  (*lib_error_handler)("String", msg);
}

String::operator const char*() const
{ 
  return (const char*)chars();
}

//  globals

StrRep  _nilStrRep = { 0, 1, { 0 } }; // nil strings point here
String _nilString;               // nil SubStrings point here




/*
 the following inline fcts are specially designed to work
 in support of String classes, and are not meant as generic replacements
 for libc "str" functions.

 inline copy fcts -  I like left-to-right from->to arguments.
 all versions assume that `to' argument is non-null

 These are worth doing inline, rather than through calls because,
 via procedural integration, adjacent copy calls can be smushed
 together by the optimizer.
*/

// copy n bytes
inline static void ncopy(const char* from, char* to, INT n)
{
  if (from != to) while (--n >= 0) *to++ = *from++;
}

// copy n bytes, null-terminate
inline static void ncopy0(const char* from, char* to, INT n)
{
  if (from != to) 
  {
    while (--n >= 0) *to++ = *from++;
    *to = 0;
  }
  else
    to[n] = 0;
}

// copy until null
inline static void scopy(const char* from, char* to)
{
  if (from != 0) while((*to++ = *from++) != 0);
}

// copy right-to-left
inline static void revcopy(const char* from, char* to, short n)
{
  if (from != 0) while (--n >= 0) *to-- = *from--;
}


inline static INT slen(const char* t) // inline  strlen
{
  if (t == 0)
    return 0;
  else
  {
    const char* a = t;
    while (*a++ != 0);
    return a - 1 - t;
  }
}

// minimum & maximum representable rep size

#define MAXStrRep_SIZE   ((1 << (sizeof(short) * CHAR_BIT - 1)) - 1)
#define MINStrRep_SIZE   16

#ifndef MALLOC_MIN_OVERHEAD
#define MALLOC_MIN_OVERHEAD  4
#endif

// The basic allocation primitive:
// Always round request to something close to a power of two.
// This ensures a bit of padding, which often means that
// concatenations don't have to realloc. Plus it tends to
// be faster when lots of Strings are created and discarded,
// since just about any version of malloc (op new()) will
// be faster when it can reuse identically-sized chunks

inline static StrRep* Snew(INT newsiz)
{
  unsigned INT siz = sizeof(StrRep) + newsiz + MALLOC_MIN_OVERHEAD;
  unsigned INT allocsiz = MINStrRep_SIZE;
  while (allocsiz < siz) allocsiz <<= 1;
  allocsiz -= MALLOC_MIN_OVERHEAD;
  if (allocsiz >= MAXStrRep_SIZE)
    (*lib_error_handler)("String", "Requested length out of range");
    
  StrRep* rep = new (operator new (allocsiz)) StrRep;
  rep->sz = allocsiz - sizeof(StrRep);
  return rep;
}

// Do-something-while-allocating routines.

// We live with two ways to signify empty Sreps: either the
// null pointer (0) or a pointer to the nilStrRep.

// We always signify unknown source lengths (usually when fed a char*)
// via len == -1, in which case it is computed.

// allocate, copying src if nonull

StrRep* Salloc(StrRep* old, const char* src, INT srclen, INT newlen)
{
  if (old == &_nilStrRep) old = 0;
  if (srclen < 0) srclen = slen(src);
  if (newlen < srclen) newlen = srclen;
  StrRep* rep;
  if (old == 0 || newlen > old->sz)
    rep = Snew(newlen);
  else
    rep = old;

  rep->len = newlen;
  ncopy0(src, rep->s, srclen);

  if (old != rep && old != 0) delete old;

  return rep;
}

// reallocate: Given the initial allocation scheme, it will
// generally be faster in the long run to get new space & copy
// than to call realloc

static StrRep*
Sresize(StrRep* old, INT newlen)
{
  if (old == &_nilStrRep) old = 0;
  StrRep* rep;
  if (old == 0)
    rep = Snew(newlen);
  else if (newlen > old->sz)
  {
    rep = Snew(newlen);
    ncopy0(old->s, rep->s, old->len);
    delete old;
  }
  else
    rep = old;

  rep->len = newlen;

  return rep;
}

void
String::alloc (INT newsize)
{
  unsigned short old_len = rep->len;
  rep = Sresize(rep, newsize);
  rep->len = old_len;
}

// like allocate, but we know that src is a StrRep

StrRep* Scopy(StrRep* old, const StrRep* s)
{
  if (old == &_nilStrRep) old = 0;
  if (s == &_nilStrRep) s = 0;
  if (old == s) 
    return (old == 0)? &_nilStrRep : old;
  else if (s == 0)
  {
    old->s[0] = 0;
    old->len = 0;
    return old;
  }
  else 
  {
    StrRep* rep;
    INT newlen = (short) s->len;	// Beware of sign extension!
    if (old == 0 || newlen > old->sz)
    {
      if (old != 0) delete old;
      rep = Snew(newlen);
    }
    else
      rep = old;
    rep->len = newlen;
    ncopy0(s->s, rep->s, newlen);
    return rep;
  }
}

// allocate & concatenate

StrRep* Scat(StrRep* old, const char* s, INT srclen, const char* t, INT tlen)
{
  if (old == &_nilStrRep) old = 0;
  if (srclen < 0) srclen = slen(s);
  if (tlen < 0) tlen = slen(t);
  INT newlen = srclen + tlen;
  StrRep* rep;

  if (old == 0 || newlen > old->sz || 
      (t >= old->s && t < &(old->s[old->len]))) // beware of aliasing
    rep = Snew(newlen);
  else
    rep = old;

  rep->len = newlen;

  ncopy(s, rep->s, srclen);
  ncopy0(t, &(rep->s[srclen]), tlen);

  if (old != rep && old != 0) delete old;

  return rep;
}

// double-concatenate

StrRep* Scat(StrRep* old, const char* s, INT srclen, const char* t, INT tlen,
             const char* u, INT ulen)
{
  if (old == &_nilStrRep) old = 0;
  if (srclen < 0) srclen = slen(s);
  if (tlen < 0) tlen = slen(t);
  if (ulen < 0) ulen = slen(u);
  INT newlen = srclen + tlen + ulen;
  StrRep* rep;
  if (old == 0 || newlen > old->sz || 
      (t >= old->s && t < &(old->s[old->len])) ||
      (u >= old->s && u < &(old->s[old->len])))
    rep = Snew(newlen);
  else
    rep = old;

  rep->len = newlen;

  ncopy(s, rep->s, srclen);
  ncopy(t, &(rep->s[srclen]), tlen);
  ncopy0(u, &(rep->s[srclen+tlen]), ulen);

  if (old != rep && old != 0) delete old;

  return rep;
}

// like cat, but we know that new stuff goes in the front of existing rep

StrRep* Sprepend(StrRep* old, const char* t, INT tlen)
{
  char* s;
  INT srclen;
  if (old == &_nilStrRep || old == 0)
  {
    s = 0; old = 0; srclen = 0;
  }
  else
  {
    s = old->s; srclen = (short) old->len;
  }
  if (tlen < 0) tlen = slen(t);
  INT newlen = srclen + tlen;
  StrRep* rep;
  if (old == 0 || newlen > old->sz || 
      (t >= old->s && t < &(old->s[old->len])))
    rep = Snew(newlen);
  else
    rep = old;

  rep->len = newlen;

  revcopy(&(s[srclen]), &(rep->s[newlen]), srclen+1);
  ncopy(t, rep->s, tlen);

  if (old != rep && old != 0) delete old;

  return rep;
}


// string compare: first argument is known to be non-null

inline static INT scmp(const char* a, const char* b)
{
  if (b == 0)
    return *a != 0;
  else
  {
    INT diff = 0;
    while ((diff = *a - *b++) == 0 && *a++ != 0);
    return diff;
  }
}


inline static INT ncmp(const char* a, INT al, const char* b, INT bl)
{
  INT n = (al <= bl)? al : bl;
  INT diff;
  while (n-- > 0) if ((diff = *a++ - *b++) != 0) return diff;
  return al - bl;
}

INT fcompare(const String& x, const String& y)
{
  const char* a = x.chars();
  const char* b = y.chars();
  INT al = x.length();
  INT bl = y.length();
  INT n = (al <= bl)? al : bl;
  INT diff = 0;
  while (n-- > 0)
  {
    char ac = *a++;
    char bc = *b++;
    if ((diff = ac - bc) != 0)
    {
      if (ac >= 'a' && ac <= 'z')
        ac = ac - 'a' + 'A';
      if (bc >= 'a' && bc <= 'z')
        bc = bc - 'a' + 'A';
      if ((diff = ac - bc) != 0)
        return diff;
    }
  }
  return al - bl;
}

// these are not inline, but pull in the above inlines, so are 
// pretty fast

INT compare(const String& x, const char* b)
{
  return scmp(x.chars(), b);
}

INT compare(const String& x, const String& y)
{
  return scmp(x.chars(), y.chars());
}

INT compare(const String& x, const SubString& y)
{
  return ncmp(x.chars(), x.length(), y.chars(), y.length());
}

INT compare(const SubString& x, const String& y)
{
  return ncmp(x.chars(), x.length(), y.chars(), y.length());
}

INT compare(const SubString& x, const SubString& y)
{
  return ncmp(x.chars(), x.length(), y.chars(), y.length());
}

INT compare(const SubString& x, const char* b)
{
  if (b == 0)
    return x.length();
  else
  {
    const char* a = x.chars();
    INT n = x.length();
    INT diff;
    while (n-- > 0) if ((diff = *a++ - *b++) != 0) return diff;
    return (*b == 0) ? 0 : -1;
  }
}

/*
 index fcts
*/

INT String::search(INT start, INT sl, char c) const
{
  const char* s = chars();
  if (sl > 0)
  {
    if (start >= 0)
    {
      const char* a = &(s[start]);
      const char* lasta = &(s[sl]);
      while (a < lasta) if (*a++ == c) return --a - s;
    }
    else
    {
      const char* a = &(s[sl + start + 1]);
      while (--a >= s) if (*a == c) return a - s;
    }
  }
  return -1;
}

INT String::search(INT start, INT sl, const char* t, INT tl) const
{
  const char* s = chars();
  if (tl < 0) tl = slen(t);
  if (sl > 0 && tl > 0)
  {
    if (start >= 0)
    {
      const char* lasts = &(s[sl - tl]);
      const char* lastt = &(t[tl]);
      const char* p = &(s[start]);

      while (p <= lasts)
      {
        const char* x = p++;
        const char* y = t;
        while (*x++ == *y++) if (y >= lastt) return --p - s;
      }
    }
    else
    {
      const char* firsts = &(s[tl - 1]);
      const char* lastt =  &(t[tl - 1]);
      const char* p = &(s[sl + start + 1]); 

      while (--p >= firsts)
      {
        const char* x = p;
        const char* y = lastt;
        while (*x-- == *y--) if (y < t) return ++x - s;
      }
    }
  }
  return -1;
}

INT String::match(INT start, INT sl, INT exact, const char* t, INT tl) const
{
  if (tl < 0) tl = slen(t);

  if (start < 0)
  {
    start = sl + start - tl + 1;
    if (start < 0 || (exact && start != 0))
      return -1;
  }
  else if (exact && sl - start != tl)
    return -1;

  if (sl == 0 || tl == 0 || sl - start < tl || start >= sl)
    return -1;

  INT n = tl;
  const char* s = &(rep->s[start]);
  while (--n >= 0) if (*s++ != *t++) return -1;
  return tl;
}

void SubString::assign(const StrRep* ysrc, const char* ys, INT ylen)
{
  if (&S == &_nilString) return;

  if (ylen < 0) ylen = slen(ys);
  StrRep* targ = S.rep;
  INT sl = (short) (targ->len) - len + ylen;

  if (ysrc == targ || sl >= targ->sz)
  {
    StrRep* oldtarg = targ;
    targ = Sresize(0, sl);
    ncopy(oldtarg->s, targ->s, pos);
    ncopy(ys, &(targ->s[pos]), ylen);
    scopy(&(oldtarg->s[pos + len]), &(targ->s[pos + ylen]));
    delete oldtarg;
  }
  else if (len == ylen)
    ncopy(ys, &(targ->s[pos]), len);
  else if (ylen < len)
  {
    ncopy(ys, &(targ->s[pos]), ylen);
    scopy(&(targ->s[pos + len]), &(targ->s[pos + ylen]));
  }
  else
  {
    revcopy(&(targ->s[targ->len]), &(targ->s[sl]), targ->len-pos-len +1);
    ncopy(ys, &(targ->s[pos]), ylen);
  }
  targ->len = sl;
  S.rep = targ;
}



/*
 * substitution
 */


INT String::_gsub(const char* pat, INT pl, const char* r, INT rl)
{
  INT nmatches = 0;
  if (pl < 0) pl = slen(pat);
  if (rl < 0) rl = slen(r);
  INT sl = length();
  if (sl <= 0 || pl <= 0 || sl < pl)
    return nmatches;
  
  const char* s = chars();

  // prepare to make new rep
  StrRep* nrep = 0;
  INT nsz = 0;
  char* x = 0;

  INT si = 0;
  INT xi = 0;
  INT remaining = sl;

  while (remaining >= pl)
  {
    INT pos = search(si, sl, pat, pl);
    if (pos < 0)
      break;
    else
    {
      ++nmatches;
      INT mustfit = xi + remaining + rl - pl;
      if (mustfit >= nsz)
      {
        if (nrep != 0) nrep->len = xi;
        nrep = Sresize(nrep, mustfit);
        nsz = (short) nrep->sz;
        x = nrep->s;
      }
      pos -= si;
      ncopy(&(s[si]), &(x[xi]), pos);
      ncopy(r, &(x[xi + pos]), rl);
      si += pos + pl;
      remaining -= pos + pl;
      xi += pos + rl;
    }
  }

  if (nrep == 0)
  {
    if (nmatches == 0)
      return nmatches;
    else
      nrep = Sresize(nrep, xi+remaining);
  }

  ncopy0(&(s[si]), &(x[xi]), remaining);
  nrep->len = xi + remaining;

  if (nrep->len <= rep->sz)   // fit back in if possible
  {
    rep->len = nrep->len;
    ncopy0(nrep->s, rep->s, rep->len);
    delete(nrep);
  }
  else
  {
    delete(rep);
    rep = nrep;
  }
  return nmatches;
}

INT String::_gsub(const Regex& pat, const char* r, INT rl)
{
  INT nmatches = 0;
  INT sl = length();
  if (sl <= 0)
    return nmatches;

  if (rl < 0) rl = slen(r);

  const char* s = chars();

  StrRep* nrep = 0;
  INT nsz = 0;

  char* x = 0;

  INT si = 0;
  INT xi = 0;
  INT remaining = sl;
  INT  pos, pl = 0;				  // how long is a regular expression?

  while (remaining > 0)
  {
    pos = pat.search(s, sl, pl, si); // unlike string search, the pos returned here is absolute
    if (pos < 0 || pl <= 0)
      break;
    else
    {
      ++nmatches;
      INT mustfit = xi + remaining + rl - pl;
      if (mustfit >= nsz)
      {
        if (nrep != 0) nrep->len = xi;
        nrep = Sresize(nrep, mustfit);
        x = nrep->s;
        nsz = (short) nrep->sz;
      }
      pos -= si;
      ncopy(&(s[si]), &(x[xi]), pos);
      ncopy(r, &(x[xi + pos]), rl);
      si += pos + pl;
      remaining -= pos + pl;
      xi += pos + rl;
    }
  }

  if (nrep == 0)
  {
    if (nmatches == 0)
      return nmatches;
    else
      nrep = Sresize(nrep, xi+remaining);
  }

  ncopy0(&(s[si]), &(x[xi]), remaining);
  nrep->len = xi + remaining;

  if (nrep->len <= rep->sz)   // fit back in if possible
  {
    rep->len = nrep->len;
    ncopy0(nrep->s, rep->s, rep->len);
    delete(nrep);
  }
  else
  {
    delete(rep);
    rep = nrep;
  }
  return nmatches;
}


/*
 * deletion
 */

void String::del(INT pos, INT len)
{
  if (pos < 0 || len <= 0 || (unsigned)(pos + len) > length()) return;
  INT nlen = length() - len;
  INT first = pos + len;
  ncopy0(&(rep->s[first]), &(rep->s[pos]), length() - first);
  rep->len = nlen;
}

void String::del(const Regex& r, INT startpos)
{
  INT mlen;
  INT first = r.search(chars(), length(), mlen, startpos);
  del(first, mlen);
}

void String::del(const char* t, INT startpos)
{
  INT tlen = slen(t);
  INT p = search(startpos, length(), t, tlen);
  del(p, tlen);
}

void String::del(const String& y, INT startpos)
{
  del(search(startpos, length(), y.chars(), y.length()), y.length());
}

void String::del(const SubString& y, INT startpos)
{
  del(search(startpos, length(), y.chars(), y.length()), y.length());
}

void String::del(char c, INT startpos)
{
  del(search(startpos, length(), c), 1);
}

/*
 * substring extraction
 */


SubString String::at(INT first, INT len)
{
  return _substr(first, len);
}

SubString String::operator() (INT first, INT len)
{
  return _substr(first, len);
}

SubString String::before(INT pos)
{
  return _substr(0, pos);
}

SubString String::through(INT pos)
{
  return _substr(0, pos+1);
}

SubString String::after(INT pos)
{
  return _substr(pos + 1, length() - (pos + 1));
}

SubString String::from(INT pos)
{
  return _substr(pos, length() - pos);
}

SubString String::at(const String& y, INT startpos)
{
  INT first = search(startpos, length(), y.chars(), y.length());
  return _substr(first,  y.length());
}

SubString String::at(const SubString& y, INT startpos)
{
  INT first = search(startpos, length(), y.chars(), y.length());
  return _substr(first, y.length());
}

SubString String::at(const Regex& r, INT startpos)
{
  INT mlen;
  INT first = r.search(chars(), length(), mlen, startpos);
  return _substr(first, mlen);
}

SubString String::at(const char* t, INT startpos)
{
  INT tlen = slen(t);
  INT first = search(startpos, length(), t, tlen);
  return _substr(first, tlen);
}

SubString String::at(char c, INT startpos)
{
  INT first = search(startpos, length(), c);
  return _substr(first, 1);
}

SubString String::before(const String& y, INT startpos)
{
  INT last = search(startpos, length(), y.chars(), y.length());
  return _substr(0, last);
}

SubString String::before(const SubString& y, INT startpos)
{
  INT last = search(startpos, length(), y.chars(), y.length());
  return _substr(0, last);
}

SubString String::before(const Regex& r, INT startpos)
{
  INT mlen;
  INT first = r.search(chars(), length(), mlen, startpos);
  return _substr(0, first);
}

SubString String::before(char c, INT startpos)
{
  INT last = search(startpos, length(), c);
  return _substr(0, last);
}

SubString String::before(const char* t, INT startpos)
{
  INT tlen = slen(t);
  INT last = search(startpos, length(), t, tlen);
  return _substr(0, last);
}

SubString String::through(const String& y, INT startpos)
{
  INT last = search(startpos, length(), y.chars(), y.length());
  if (last >= 0) last += y.length();
  return _substr(0, last);
}

SubString String::through(const SubString& y, INT startpos)
{
  INT last = search(startpos, length(), y.chars(), y.length());
  if (last >= 0) last += y.length();
  return _substr(0, last);
}

SubString String::through(const Regex& r, INT startpos)
{
  INT mlen;
  INT first = r.search(chars(), length(), mlen, startpos);
  if (first >= 0) first += mlen;
  return _substr(0, first);
}

SubString String::through(char c, INT startpos)
{
  INT last = search(startpos, length(), c);
  if (last >= 0) last += 1;
  return _substr(0, last);
}

SubString String::through(const char* t, INT startpos)
{
  INT tlen = slen(t);
  INT last = search(startpos, length(), t, tlen);
  if (last >= 0) last += tlen;
  return _substr(0, last);
}

SubString String::after(const String& y, INT startpos)
{
  INT first = search(startpos, length(), y.chars(), y.length());
  if (first >= 0) first += y.length();
  return _substr(first, length() - first);
}

SubString String::after(const SubString& y, INT startpos)
{
  INT first = search(startpos, length(), y.chars(), y.length());
  if (first >= 0) first += y.length();
  return _substr(first, length() - first);
}

SubString String::after(char c, INT startpos)
{
  INT first = search(startpos, length(), c);
  if (first >= 0) first += 1;
  return _substr(first, length() - first);
}

SubString String::after(const Regex& r, INT startpos)
{
  INT mlen;
  INT first = r.search(chars(), length(), mlen, startpos);
  if (first >= 0) first += mlen;
  return _substr(first, length() - first);
}

SubString String::after(const char* t, INT startpos)
{
  INT tlen = slen(t);
  INT first = search(startpos, length(), t, tlen);
  if (first >= 0) first += tlen;
  return _substr(first, length() - first);
}

SubString String::from(const String& y, INT startpos)
{
  INT first = search(startpos, length(), y.chars(), y.length());
  return _substr(first, length() - first);
}

SubString String::from(const SubString& y, INT startpos)
{
  INT first = search(startpos, length(), y.chars(), y.length());
  return _substr(first, length() - first);
}

SubString String::from(const Regex& r, INT startpos)
{
  INT mlen;
  INT first = r.search(chars(), length(), mlen, startpos);
  return _substr(first, length() - first);
}

SubString String::from(char c, INT startpos)
{
  INT first = search(startpos, length(), c);
  return _substr(first, length() - first);
}

SubString String::from(const char* t, INT startpos)
{
  INT tlen = slen(t);
  INT first = search(startpos, length(), t, tlen);
  return _substr(first, length() - first);
}



/*
 * split/join
 */


INT split(const String& src, String results[], INT n, const String& sep)
{
  String x = src;
  const char* s = x.chars();
  INT sl = x.length();
  INT i = 0;
  INT pos = 0;
  while (i < n && pos < sl)
  {
    INT p = x.search(pos, sl, sep.chars(), sep.length());
    if (p < 0)
      p = sl;
    results[i].rep = Salloc(results[i].rep, &(s[pos]), p - pos, p - pos);
    i++;
    pos = p + sep.length();
  }
  return i;
}

INT split(const String& src, String results[], INT n, const Regex& r)
{
  String x = src;
  const char* s = x.chars();
  INT sl = x.length();
  INT i = 0;
  INT pos = 0;
  INT p, matchlen;
  while (i < n && pos < sl)
  {
    p = r.search(s, sl, matchlen, pos);
    if (p < 0)
      p = sl;
    results[i].rep = Salloc(results[i].rep, &(s[pos]), p - pos, p - pos);
    i++;
    pos = p + matchlen;
  }
  return i;
}


#if defined(__GNUG__) && !defined(_G_NO_NRV)
#define RETURN(r) return
#define RETURNS(r) return r;
#define RETURN_OBJECT(TYPE, NAME) /* nothing */
#else /* _G_NO_NRV */
#define RETURN(r) return r
#define RETURNS(r) /* nothing */
#define RETURN_OBJECT(TYPE, NAME) TYPE NAME;
#endif

String join(String src[], INT n, const String& separator) RETURNS(x)
{
  RETURN_OBJECT(String,x)
  String sep = separator;
  INT xlen = 0;
  INT i;
  for (i = 0; i < n; ++i)
    xlen += src[i].length();
  xlen += (n - 1) * sep.length();

  x.rep = Sresize (x.rep, xlen);

  INT j = 0;
  
  for (i = 0; i < n - 1; ++i)
  {
    ncopy(src[i].chars(), &(x.rep->s[j]), src[i].length());
    j += src[i].length();
    ncopy(sep.chars(), &(x.rep->s[j]), sep.length());
    j += sep.length();
  }
  ncopy0(src[i].chars(), &(x.rep->s[j]), src[i].length());
  RETURN(x);
}
  
/*
 misc
*/

    
StrRep* Sreverse(const StrRep* src, StrRep* dest)
{
  INT n = (short) src->len;
  if (src != dest)
    dest = Salloc(dest, src->s, n, n);
  if (n > 0)
  {
    char* a = dest->s;
    char* b = &(a[n - 1]);
    while (a < b)
    {
      char t = *a;
      *a++ = *b;
      *b-- = t;
    }
  }
  return dest;
}


StrRep* Supcase(const StrRep* src, StrRep* dest)
{
  INT n = (short) src->len;
  if (src != dest) dest = Salloc(dest, src->s, n, n);
  char* p = dest->s;
  char* e = &(p[n]);
  for (; p < e; ++p) if (islower(*p)) *p = toupper(*p);
  return dest;
}

StrRep* Sdowncase(const StrRep* src, StrRep* dest)
{
  INT n = (short) src->len;
  if (src != dest) dest = Salloc(dest, src->s, n, n);
  char* p = dest->s;
  char* e = &(p[n]);
  for (; p < e; ++p) if (isupper(*p)) *p = tolower(*p);
  return dest;
}

StrRep* Scapitalize(const StrRep* src, StrRep* dest)
{
  INT n = (short) src->len;
  if (src != dest) dest = Salloc(dest, src->s, n, n);

  char* p = dest->s;
  char* e = &(p[n]);
  for (; p < e; ++p)
  {
    INT at_word = islower(*p);

    if (at_word)
      *p = toupper(*p);
    else 
      at_word = isupper(*p) || isdigit(*p);

    if (at_word)
    {
      while (++p < e)
      {
        if (isupper(*p))
          *p = tolower(*p);
	/* A '\'' does not break a word, so that "Nathan's" stays
	   "Nathan's" rather than turning into "Nathan'S". */
        else if (!islower(*p) && !isdigit(*p) && (*p != '\''))
          break;
      }
    }
  }
  return dest;
}

#if defined(__GNUG__) && !defined(_G_NO_NRV)

String replicate(char c, INT n) return w;
{
  w.rep = Sresize(w.rep, n);
  char* p = w.rep->s;
  while (n-- > 0) *p++ = c;
  *p = 0;
}

String replicate(const String& y, INT n) return w
{
  INT len = y.length();
  w.rep = Sresize(w.rep, n * len);
  char* p = w.rep->s;
  while (n-- > 0)
  {
    ncopy(y.chars(), p, len);
    p += len;
  }
  *p = 0;
}

String common_prefix(const String& x, const String& y, INT startpos) return r;
{
  const char* xchars = x.chars();
  const char* ychars = y.chars();
  const char* xs = &(xchars[startpos]);
  const char* ss = xs;
  const char* topx = &(xchars[x.length()]);
  const char* ys = &(ychars[startpos]);
  const char* topy = &(ychars[y.length()]);
  INT l;
  for (l = 0; xs < topx && ys < topy && *xs++ == *ys++; ++l);
  r.rep = Salloc(r.rep, ss, l, l);
}

String common_suffix(const String& x, const String& y, INT startpos) return r;
{
  const char* xchars = x.chars();
  const char* ychars = y.chars();
  const char* xs = &(xchars[x.length() + startpos]);
  const char* botx = xchars;
  const char* ys = &(ychars[y.length() + startpos]);
  const char* boty = ychars;
  INT l;
  for (l = 0; xs >= botx && ys >= boty && *xs == *ys ; --xs, --ys, ++l);
  r.rep = Salloc(r.rep, ++xs, l, l);
}

#else

String replicate(char c, INT n)
{
  String w;
  w.rep = Sresize(w.rep, n);
  char* p = w.rep->s;
  while (n-- > 0) *p++ = c;
  *p = 0;
  return w;
}

String replicate(const String& y, INT n)
{
  String w;
  INT len = y.length();
  w.rep = Sresize(w.rep, n * len);
  char* p = w.rep->s;
  while (n-- > 0)
  {
    ncopy(y.chars(), p, len);
    p += len;
  }
  *p = 0;
  return w;
}

String common_prefix(const String& x, const String& y, INT startpos)
{
  String r;
  const char* xchars = x.chars();
  const char* ychars = y.chars();
  const char* xs = &(xchars[startpos]);
  const char* ss = xs;
  const char* topx = &(xchars[x.length()]);
  const char* ys = &(ychars[startpos]);
  const char* topy = &(ychars[y.length()]);
  INT l;
  for (l = 0; xs < topx && ys < topy && *xs++ == *ys++; ++l);
  r.rep = Salloc(r.rep, ss, l, l);
  return r;
}

String common_suffix(const String& x, const String& y, INT startpos) 
{
  String r;
  const char* xchars = x.chars();
  const char* ychars = y.chars();
  const char* xs = &(xchars[x.length() + startpos]);
  const char* botx = xchars;
  const char* ys = &(ychars[y.length() + startpos]);
  const char* boty = ychars;
  INT l;
  for (l = 0; xs >= botx && ys >= boty && *xs == *ys ; --xs, --ys, ++l);
  r.rep = Salloc(r.rep, ++xs, l, l);
  return r;
}

#endif

// IO

istream& operator>>(istream& s, String& x)
{

  istream::sentry cerberos(s);
  if (!cerberos || (!(s.flags() & ios::skipws) && !ws(s)))
  {
    s.clear(ios::failbit|s.rdstate()); // Redundant if using GNU iostreams.
    return s;
  }
  INT ch;
  INT i = 0;
  x.rep = Sresize(x.rep, 20);
  register streambuf *sb = s.rdbuf();
  while ((ch = sb->sbumpc()) != EOF)
  {
    if (isspace(ch))
      break;
    if (i >= x.rep->sz - 1)
      x.rep = Sresize(x.rep, i+1);
    x.rep->s[i++] = ch;
  }
  x.rep->s[i] = 0;
  x.rep->len = i;
  ios_base::iostate new_state = s.rdstate();
  if (i == 0) new_state |= ios::failbit;
  if (ch == EOF) new_state |= ios::eofbit;
  s.clear(new_state);
  return s;
}

INT readline(istream& s, String& x, char terminator, INT discard)
{
  istream::sentry cerberos(s);
  if (!cerberos)
    return 0;
  INT ch;
  INT i = 0;
  x.rep = Sresize(x.rep, 80);
  register streambuf *sb = s.rdbuf();
  while ((ch = sb->sbumpc()) != EOF)
  {
    if (ch != terminator || !discard)
    {
      if (i >= x.rep->sz - 1)
        x.rep = Sresize(x.rep, i+1);
      x.rep->s[i++] = ch;
    }
    if (ch == terminator)
      break;
  }
  x.rep->s[i] = 0;
  x.rep->len = i;
  if (ch == EOF) s.clear(ios::eofbit|s.rdstate());
  return i;
}


ostream& operator<<(ostream& s, const SubString& x)
{ 
  const char* a = x.chars();
  const char* lasta = &(a[x.length()]);
  while (a < lasta)
    s.put(*a++);
  return(s);
}

// from John.Willis@FAS.RI.CMU.EDU

INT String::freq(const SubString& y) const
{
  INT found = 0;
  for (unsigned INT i = 0; i < length(); i++) 
    if (match(i,length(),0,y.chars(), y.length())>= 0) found++;
  return(found);
}

INT String::freq(const String& y) const
{
  INT found = 0;
  for (unsigned INT i = 0; i < length(); i++) 
    if (match(i,length(),0,y.chars(),y.length()) >= 0) found++;
  return(found);
}

INT String::freq(const char* t) const
{
  INT found = 0;
  for (unsigned INT i = 0; i < length(); i++) 
    if (match(i,length(),0,t) >= 0) found++;
  return(found);
}

INT String::freq(char c) const
{
  INT found = 0;
  for (unsigned INT i = 0; i < length(); i++) 
    if (match(i,length(),0,&c,1) >= 0) found++;
  return(found);
}


INT String::OK() const
{
  if (rep == 0             // don't have a rep
    || rep->len > rep->sz     // string oustide bounds
    || rep->s[rep->len] != 0)   // not null-terminated
      error("invariant failure");
  return 1;
}

INT SubString::OK() const
{
  INT v = S != (const char*)0; // have a String;
  v &= S.OK();                 // that is legal
  v &= pos + len >= (short) S.rep->len;// pos and len within bounds
  if (!v) S.error("SubString invariant failure");
  return v;
}

