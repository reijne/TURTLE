// This may look like C code, but it is really -*- C++ -*-
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

#ifdef __GNUG__
#define _G_NO_NRV
#endif

#ifndef _String_h
#define _String_h 1

#include <iostream>

#include "../../config.h"

#include "Regex.h"

#undef OK

struct StrRep                     // internal String representations
{
  unsigned short    len;         // string length 
  unsigned short    sz;          // allocated space
  char              s[1];        // the string starts here 
                                 // (at least 1 char for trailing null)
                                 // allocated & expanded via non-public fcts
};

// primitive ops on StrReps -- nearly all String fns go through these.

StrRep*     Salloc(StrRep*, const char*, INT, INT);
StrRep*     Scopy(StrRep*, const StrRep*);
StrRep*     Scat(StrRep*, const char*, INT, const char*, INT);
StrRep*     Scat(StrRep*, const char*, INT,const char*,INT, const char*,INT);
StrRep*     Sprepend(StrRep*, const char*, INT);
StrRep*     Sreverse(const StrRep*, StrRep*);
StrRep*     Supcase(const StrRep*, StrRep*);
StrRep*     Sdowncase(const StrRep*, StrRep*);
StrRep*     Scapitalize(const StrRep*, StrRep*);

// These classes need to be defined in the order given

class String;
class SubString;

class SubString
{

  friend class      String;
protected:

  String&           S;        // The String I'm a substring of
  unsigned short    pos;      // starting position in S's rep
  unsigned short    len;      // length of substring

  void              assign(const StrRep*, const char*, INT = -1);

public:
                    SubString(String& x, INT p, INT l);
		    SubString(const SubString& x);
                   ~SubString();

  SubString&        operator =  (const String&     y);
  SubString&        operator =  (const SubString&  y);
  SubString&        operator =  (const char* t);
  SubString&        operator =  (char        c);

// return 1 if target appears anywhere in SubString; else 0

  INT               contains(char        c) const;
  INT               contains(const String&     y) const;
  INT               contains(const SubString&  y) const;
  INT               contains(const char* t) const;
  INT               contains(const Regex&       r) const;

// return 1 if target matches entire SubString

  INT               matches(const Regex&  r) const;

// IO 
  friend std::ostream&   operator<<(std::ostream & s, const SubString& x);

// status

  unsigned INT      length() const;
  INT               empty() const;
  const char*       chars() const;

  INT               OK() const; 

};


class String
{
  friend class      SubString;

protected:
  StrRep*           rep;   // Strings are pointers to their representations

// some helper functions

  INT               search(INT, INT, const char*, INT = -1) const;
  INT               search(INT, INT, char) const;
  INT               match(INT, INT, INT, const char*, INT = -1) const;
  INT               _gsub(const char*, INT, const char* ,INT);
  INT               _gsub(const Regex&, const char*, INT);
  SubString         _substr(INT, INT);

public:

// constructors & assignment

                    String();
                    String(const String& x);
                    String(const SubString&  x);
                    String(const char* t);
                    String(const char* t, INT len);
                    String(char c);

                    ~String();

  String&           operator =  (const String&     y);
  String&           operator =  (const char* y);
  String&           operator =  (char        c);
  String&           operator =  (const SubString&  y);

// concatenation

  String&           operator += (const String&     y); 
  String&           operator += (const SubString&  y);
  String&           operator += (const char* t);
  String&           operator += (char        c);

  void              prepend(const String&     y); 
  void              prepend(const SubString&  y);
  void              prepend(const char* t);
  void              prepend(char        c);


// procedural versions:
// concatenate first 2 args, store result in last arg

  friend inline void     cat(const String&, const String&, String&);
  friend inline void     cat(const String&, const SubString&, String&);
  friend inline void     cat(const String&, const char*, String&);
  friend inline void     cat(const String&, char, String&);

  friend inline void     cat(const SubString&, const String&, String&);
  friend inline void     cat(const SubString&, const SubString&, String&);
  friend inline void     cat(const SubString&, const char*, String&);
  friend inline void     cat(const SubString&, char, String&);

  friend inline void     cat(const char*, const String&, String&);
  friend inline void     cat(const char*, const SubString&, String&);
  friend inline void     cat(const char*, const char*, String&);
  friend inline void     cat(const char*, char, String&);

// double concatenation, by request. (yes, there are too many versions, 
// but if one is supported, then the others should be too...)
// Concatenate first 3 args, store in last arg

  friend inline void     cat(const String&,const String&, const String&,String&);
  friend inline void     cat(const String&,const String&,const SubString&,String&);
  friend inline void     cat(const String&,const String&, const char*, String&);
  friend inline void     cat(const String&,const String&, char, String&);
  friend inline void     cat(const String&,const SubString&,const String&,String&);
  inline friend void     cat(const String&,const SubString&,const SubString&,String&);
  friend inline void     cat(const String&,const SubString&, const char*, String&);
  friend inline void     cat(const String&,const SubString&, char, String&);
  friend inline void     cat(const String&,const char*, const String&,    String&);
  friend inline void     cat(const String&,const char*, const SubString&, String&);
  friend inline void     cat(const String&,const char*, const char*, String&);
  friend inline void     cat(const String&,const char*, char, String&);

  friend inline void     cat(const char*, const String&, const String&,String&);
  friend inline void     cat(const char*,const String&,const SubString&,String&);
  friend inline void     cat(const char*,const String&, const char*, String&);
  friend inline void     cat(const char*,const String&, char, String&);
  friend inline void     cat(const char*,const SubString&,const String&,String&);
  friend inline void     cat(const char*,const SubString&,const SubString&,String&);
  friend inline void     cat(const char*,const SubString&, const char*, String&);
  friend inline void     cat(const char*,const SubString&, char, String&);
  friend inline void     cat(const char*,const char*, const String&,    String&);
  friend inline void     cat(const char*,const char*, const SubString&, String&);
  friend inline void     cat(const char*,const char*, const char*, String&);
  friend inline void     cat(const char*,const char*, char, String&);


// searching & matching

// return position of target in string or -1 for failure

  INT               index(char        c, INT startpos = 0) const;      
  INT               index(const String&     y, INT startpos = 0) const;      
  INT               index(const SubString&  y, INT startpos = 0) const;      
  INT               index(const char* t, INT startpos = 0) const;  
  INT               index(const Regex&      r, INT startpos = 0) const;       

// return 1 if target appears anyhere in String; else 0

  INT               contains(char        c) const;
  INT               contains(const String&     y) const;
  INT               contains(const SubString&  y) const;
  INT               contains(const char* t) const;
  INT               contains(const Regex&      r) const;

// return 1 if target appears anywhere after position pos 
// (or before, if pos is negative) in String; else 0

  INT               contains(char        c, INT pos) const;
  INT               contains(const String&     y, INT pos) const;
  INT               contains(const SubString&  y, INT pos) const;
  INT               contains(const char* t, INT pos) const;
  INT               contains(const Regex&      r, INT pos) const;

// return 1 if target appears at position pos in String; else 0

  INT               matches(char        c, INT pos = 0) const;
  INT               matches(const String&     y, INT pos = 0) const;
  INT               matches(const SubString&  y, INT pos = 0) const;
  INT               matches(const char* t, INT pos = 0) const;
  INT               matches(const Regex&      r, INT pos = 0) const;

//  return number of occurences of target in String

  INT               freq(char        c) const; 
  INT               freq(const String&     y) const;
  INT               freq(const SubString&  y) const;
  INT               freq(const char* t) const;

// SubString extraction

// Note that you can't take a substring of a const String, since
// this leaves open the possiblility of indirectly modifying the
// String through the SubString

  SubString         at(INT         pos, INT len);
  SubString         operator () (INT         pos, INT len); // synonym for at

  SubString         at(const String&     x, INT startpos = 0); 
  SubString         at(const SubString&  x, INT startpos = 0); 
  SubString         at(const char* t, INT startpos = 0);
  SubString         at(char        c, INT startpos = 0);
  SubString         at(const Regex&      r, INT startpos = 0); 

  SubString         before(INT          pos);
  SubString         before(const String&      x, INT startpos = 0);
  SubString         before(const SubString&   x, INT startpos = 0);
  SubString         before(const char*  t, INT startpos = 0);
  SubString         before(char         c, INT startpos = 0);
  SubString         before(const Regex&       r, INT startpos = 0);

  SubString         through(INT          pos);
  SubString         through(const String&      x, INT startpos = 0);
  SubString         through(const SubString&   x, INT startpos = 0);
  SubString         through(const char*  t, INT startpos = 0);
  SubString         through(char         c, INT startpos = 0);
  SubString         through(const Regex&       r, INT startpos = 0);

  SubString         from(INT          pos);
  SubString         from(const String&      x, INT startpos = 0);
  SubString         from(const SubString&   x, INT startpos = 0);
  SubString         from(const char*  t, INT startpos = 0);
  SubString         from(char         c, INT startpos = 0);
  SubString         from(const Regex&       r, INT startpos = 0);

  SubString         after(INT         pos);
  SubString         after(const String&     x, INT startpos = 0);
  SubString         after(const SubString&  x, INT startpos = 0);
  SubString         after(const char* t, INT startpos = 0);
  SubString         after(char        c, INT startpos = 0);
  SubString         after(const Regex&      r, INT startpos = 0);


// deletion

// delete len chars starting at pos
  void              del(INT         pos, INT len);

// delete the first occurrence of target after startpos

  void              del(const String&     y, INT startpos = 0);
  void              del(const SubString&  y, INT startpos = 0);
  void              del(const char* t, INT startpos = 0);
  void              del(char        c, INT startpos = 0);
  void              del(const Regex&      r, INT startpos = 0);

// global substitution: substitute all occurrences of pat with repl

  INT               gsub(const String&     pat, const String&     repl);
  INT               gsub(const SubString&  pat, const String&     repl);
  INT               gsub(const char* pat, const String&     repl);
  INT               gsub(const char* pat, const char* repl);
  INT               gsub(const Regex&      pat, const String&     repl);

// friends & utilities

// split string into array res at separators; return number of elements

  friend INT        split(const String& x, String res[], INT maxn, 
                          const String& sep);
  friend INT        split(const String& x, String res[], INT maxn, 
                          const Regex&  sep);

  friend String     common_prefix(const String& x, const String& y, 
                                  INT startpos = 0);
  friend String     common_suffix(const String& x, const String& y, 
                                  INT startpos = -1);
  friend String     replicate(char        c, INT n);
  friend String     replicate(const String&     y, INT n);
  friend String     join(String src[], INT n, const String& sep);

// simple builtin transformations

  friend inline String     reverse(const String& x);
  friend inline String     upcase(const String& x);
  friend inline String     downcase(const String& x);
  friend inline String     capitalize(const String& x);

// in-place versions of above

  void              reverse();
  void              upcase();
  void              downcase();
  void              capitalize();

// element extraction

  char&             operator [] (INT i);
  const char&       operator [] (INT i) const;
  char              elem(INT i) const;
  char              firstchar() const;
  char              lastchar() const;

// conversion

                    operator const char*() const;
  const char*       chars() const;


// IO
  friend inline std::ostream&   operator<<(std::ostream& s, const String& x);
  friend std::ostream&   operator<<(std::ostream& s, const SubString& x);
  friend std::istream&   operator>>(std::istream& s, String& x);

  friend INT        readline(std::istream& s, String& x, 
                             char terminator = '\n',
                             INT discard_terminator = 1);
// status

  unsigned INT      length() const;
  INT               empty() const;

// preallocate some space for String
  void              alloc(INT newsize);

// report current allocation (not length!)

  INT               allocation() const;


  void     error(const char* msg) const;

  INT               OK() const;
};

typedef String StrTmp; // for backward compatibility

// other externs

INT        compare(const String&    x, const String&     y);
INT        compare(const String&    x, const SubString&  y);
INT        compare(const String&    x, const char* y);
INT        compare(const SubString& x, const String&     y);
INT        compare(const SubString& x, const SubString&  y);
INT        compare(const SubString& x, const char* y);
INT        fcompare(const String&   x, const String&     y); // ignore case

extern StrRep  _nilStrRep;
extern String _nilString;

// status reports, needed before defining other things

inline unsigned INT String::length() const {  return rep->len; }
inline INT         String::empty() const { return rep->len == 0; }
inline const char* String::chars() const { return &(rep->s[0]); }
inline INT         String::allocation() const { return rep->sz; }

inline unsigned INT SubString::length() const { return len; }
inline INT         SubString::empty() const { return len == 0; }
inline const char* SubString::chars() const { return &(S.rep->s[pos]); }


// constructors

inline String::String() 
  : rep(&_nilStrRep) {}
inline String::String(const String& x) 
  : rep(Scopy(0, x.rep)) {}
inline String::String(const char* t) 
  : rep(Salloc(0, t, -1, -1)) {}
inline String::String(const char* t, INT tlen)
  : rep(Salloc(0, t, tlen, tlen)) {}
inline String::String(const SubString& y)
  : rep(Salloc(0, y.chars(), y.length(), y.length())) {}
inline String::String(char c) 
  : rep(Salloc(0, &c, 1, 1)) {}

inline String::~String() { if (rep != &_nilStrRep) delete rep; }

inline SubString::SubString(const SubString& x)
  :S(x.S), pos(x.pos), len(x.len) {}
inline SubString::SubString(String& x, INT first, INT l)
  :S(x), pos(first), len(l) {}

inline SubString::~SubString() {}

// assignment

inline String& String::operator =  (const String& y)
{ 
  rep = Scopy(rep, y.rep);
  return *this;
}

inline String& String::operator=(const char* t)
{
  rep = Salloc(rep, t, -1, -1);
  return *this;
}

inline String& String::operator=(const SubString&  y)
{
  rep = Salloc(rep, y.chars(), y.length(), y.length());
  return *this;
}

inline String& String::operator=(char c)
{
  rep = Salloc(rep, &c, 1, 1);
  return *this;
}


inline SubString& SubString::operator = (const char* ys)
{
  assign(0, ys);
  return *this;
}

inline SubString& SubString::operator = (char ch)
{
  assign(0, &ch, 1);
  return *this;
}

inline SubString& SubString::operator = (const String& y)
{
  assign(y.rep, y.chars(), y.length());
  return *this;
}

inline SubString& SubString::operator = (const SubString& y)
{
  assign(y.S.rep, y.chars(), y.length());
  return *this;
}

// Zillions of cats...

inline void cat(const String& x, const String& y, String& r)
{
  r.rep = Scat(r.rep, x.chars(), x.length(), y.chars(), y.length());
}

inline void cat(const String& x, const SubString& y, String& r)
{
  r.rep = Scat(r.rep, x.chars(), x.length(), y.chars(), y.length());
}

inline void cat(const String& x, const char* y, String& r)
{
  r.rep = Scat(r.rep, x.chars(), x.length(), y, -1);
}

inline void cat(const String& x, char y, String& r)
{
  r.rep = Scat(r.rep, x.chars(), x.length(), &y, 1);
}

inline void cat(const SubString& x, const String& y, String& r)
{
  r.rep = Scat(r.rep, x.chars(), x.length(), y.chars(), y.length());
}

inline void cat(const SubString& x, const SubString& y, String& r)
{
  r.rep = Scat(r.rep, x.chars(), x.length(), y.chars(), y.length());
}

inline void cat(const SubString& x, const char* y, String& r)
{
  r.rep = Scat(r.rep, x.chars(), x.length(), y, -1);
}

inline void cat(const SubString& x, char y, String& r)
{
  r.rep = Scat(r.rep, x.chars(), x.length(), &y, 1);
}

inline void cat(const char* x, const String& y, String& r)
{
  r.rep = Scat(r.rep, x, -1, y.chars(), y.length());
}

inline void cat(const char* x, const SubString& y, String& r)
{
  r.rep = Scat(r.rep, x, -1, y.chars(), y.length());
}

inline void cat(const char* x, const char* y, String& r)
{
  r.rep = Scat(r.rep, x, -1, y, -1);
}

inline void cat(const char* x, char y, String& r)
{
  r.rep = Scat(r.rep, x, -1, &y, 1);
}

inline void cat(const String& a, const String& x, const String& y, String& r)
{
  r.rep = Scat(r.rep, a.chars(), a.length(), x.chars(), x.length(), y.chars(), y.length());
}

inline void cat(const String& a, const String& x, const SubString& y, String& r)
{
  r.rep = Scat(r.rep, a.chars(), a.length(), x.chars(), x.length(), y.chars(), y.length());
}

inline void cat(const String& a, const String& x, const char* y, String& r)
{
  r.rep = Scat(r.rep, a.chars(), a.length(), x.chars(), x.length(), y, -1);
}

inline void cat(const String& a, const String& x, char y, String& r)
{
  r.rep = Scat(r.rep, a.chars(), a.length(), x.chars(), x.length(), &y, 1);
}

inline void cat(const String& a, const SubString& x, const String& y, String& r)
{
  r.rep = Scat(r.rep, a.chars(), a.length(), x.chars(), x.length(), y.chars(), y.length());
}

inline void cat(const String& a, const SubString& x, const SubString& y, String& r)
{
  r.rep = Scat(r.rep, a.chars(), a.length(), x.chars(), x.length(), y.chars(), y.length());
}

inline void cat(const String& a, const SubString& x, const char* y, String& r)
{
  r.rep = Scat(r.rep, a.chars(), a.length(), x.chars(), x.length(), y, -1);
}

inline void cat(const String& a, const SubString& x, char y, String& r)
{
  r.rep = Scat(r.rep, a.chars(), a.length(), x.chars(), x.length(), &y, 1);
}

inline void cat(const String& a, const char* x, const String& y, String& r)
{
  r.rep = Scat(r.rep, a.chars(), a.length(), x, -1, y.chars(), y.length());
}

inline void cat(const String& a, const char* x, const SubString& y, String& r)
{
  r.rep = Scat(r.rep, a.chars(), a.length(), x, -1, y.chars(), y.length());
}

inline void cat(const String& a, const char* x, const char* y, String& r)
{
  r.rep = Scat(r.rep, a.chars(), a.length(), x, -1, y, -1);
}

inline void cat(const String& a, const char* x, char y, String& r)
{
  r.rep = Scat(r.rep, a.chars(), a.length(), x, -1, &y, 1);
}


inline void cat(const char* a, const String& x, const String& y, String& r)
{
  r.rep = Scat(r.rep, a, -1, x.chars(), x.length(), y.chars(), y.length());
}

inline void cat(const char* a, const String& x, const SubString& y, String& r)
{
  r.rep = Scat(r.rep, a, -1, x.chars(), x.length(), y.chars(), y.length());
}

inline void cat(const char* a, const String& x, const char* y, String& r)
{
  r.rep = Scat(r.rep, a, -1, x.chars(), x.length(), y, -1);
}

inline void cat(const char* a, const String& x, char y, String& r)
{
  r.rep = Scat(r.rep, a, -1, x.chars(), x.length(), &y, 1);
}

inline void cat(const char* a, const SubString& x, const String& y, String& r)
{
  r.rep = Scat(r.rep, a, -1, x.chars(), x.length(), y.chars(), y.length());
}

inline void cat(const char* a, const SubString& x, const SubString& y, String& r)
{
  r.rep = Scat(r.rep, a, -1, x.chars(), x.length(), y.chars(), y.length());
}

inline void cat(const char* a, const SubString& x, const char* y, String& r)
{
  r.rep = Scat(r.rep, a, -1, x.chars(), x.length(), y, -1);
}

inline void cat(const char* a, const SubString& x, char y, String& r)
{
  r.rep = Scat(r.rep, a, -1, x.chars(), x.length(), &y, 1);
}

inline void cat(const char* a, const char* x, const String& y, String& r)
{
  r.rep = Scat(r.rep, a, -1, x, -1, y.chars(), y.length());
}

inline void cat(const char* a, const char* x, const SubString& y, String& r)
{
  r.rep = Scat(r.rep, a, -1, x, -1, y.chars(), y.length());
}

inline void cat(const char* a, const char* x, const char* y, String& r)
{
  r.rep = Scat(r.rep, a, -1, x, -1, y, -1);
}

inline void cat(const char* a, const char* x, char y, String& r)
{
  r.rep = Scat(r.rep, a, -1, x, -1, &y, 1);
}


// operator versions

inline String& String::operator +=(const String& y)
{
  cat(*this, y, *this);
  return *this;
}

inline String& String::operator +=(const SubString& y)
{
  cat(*this, y, *this);
  return *this;
}

inline String& String::operator += (const char* y)
{
  cat(*this, y, *this);
  return *this;
}

inline String& String:: operator +=(char y)
{
  cat(*this, y, *this);
  return *this;
}

// constructive concatenation

#if defined(__GNUG__) && !defined(_G_NO_NRV)

inline String operator + (const String& x, const String& y) return r;
{
  cat(x, y, r);
}

inline String operator + (const String& x, const SubString& y) return r;
{
  cat(x, y, r);
}

inline String operator + (const String& x, const char* y) return r;
{
  cat(x, y, r);
}

inline String operator + (const String& x, char y) return r;
{
  cat(x, y, r);
}

inline String operator + (const SubString& x, const String& y) return r;
{
  cat(x, y, r);
}

inline String operator + (const SubString& x, const SubString& y) return r;
{
  cat(x, y, r);
}

inline String operator + (const SubString& x, const char* y) return r;
{
  cat(x, y, r);
}

inline String operator + (const SubString& x, char y) return r;
{
  cat(x, y, r);
}

inline String operator + (const char* x, const String& y) return r;
{
  cat(x, y, r);
}

inline String operator + (const char* x, const SubString& y) return r;
{
  cat(x, y, r);
}

inline String reverse(const String& x) return r;
{
  r.rep = Sreverse(x.rep, r.rep);
}

inline String upcase(const String& x) return r;
{
  r.rep = Supcase(x.rep, r.rep);
}

inline String downcase(const String& x) return r;
{
  r.rep = Sdowncase(x.rep, r.rep);
}

inline String capitalize(const String& x) return r;
{
  r.rep = Scapitalize(x.rep, r.rep);
}

#else /* NO_NRV */

inline String operator + (const String& x, const String& y)
{
  String r;  cat(x, y, r);  return r;
}

inline String operator + (const String& x, const SubString& y) 
{
  String r; cat(x, y, r); return r;
}

inline String operator + (const String& x, const char* y) 
{
  String r; cat(x, y, r); return r;
}

inline String operator + (const String& x, char y) 
{
  String r; cat(x, y, r); return r;
}

inline String operator + (const SubString& x, const String& y) 
{
  String r; cat(x, y, r); return r;
}

inline String operator + (const SubString& x, const SubString& y) 
{
  String r; cat(x, y, r); return r;
}

inline String operator + (const SubString& x, const char* y) 
{
  String r; cat(x, y, r); return r;
}

inline String operator + (const SubString& x, char y) 
{
  String r; cat(x, y, r); return r;
}

inline String operator + (const char* x, const String& y) 
{
  String r; cat(x, y, r); return r;
}

inline String operator + (const char* x, const SubString& y) 
{
  String r; cat(x, y, r); return r;
}

inline String reverse(const String& x) 
{
  String r; r.rep = Sreverse(x.rep, r.rep); return r;
}

inline String upcase(const String& x) 
{
  String r; r.rep = Supcase(x.rep, r.rep); return r;
}

inline String downcase(const String& x) 
{
  String r; r.rep = Sdowncase(x.rep, r.rep); return r;
}

inline String capitalize(const String& x) 
{
  String r; r.rep = Scapitalize(x.rep, r.rep); return r;
}

#endif

// prepend

inline void String::prepend(const String& y)
{
  rep = Sprepend(rep, y.chars(), y.length());
}

inline void String::prepend(const char* y)
{
  rep = Sprepend(rep, y, -1); 
}

inline void String::prepend(char y)
{
  rep = Sprepend(rep, &y, 1); 
}

inline void String::prepend(const SubString& y)
{
  rep = Sprepend(rep, y.chars(), y.length());
}

// misc transformations


inline void String::reverse()
{
  rep = Sreverse(rep, rep);
}


inline void String::upcase()
{
  rep = Supcase(rep, rep);
}


inline void String::downcase()
{
  rep = Sdowncase(rep, rep);
}


inline void String::capitalize()
{
  rep = Scapitalize(rep, rep);
}

// element extraction

inline char&  String::operator [] (INT i) 
{ 
  if (((unsigned)i) >= length()) error("invalid index");
  return rep->s[i];
}

inline const char&  String::operator [] (INT i) const
{ 
  if (((unsigned)i) >= length()) error("invalid index");
  return rep->s[i];
}

inline char  String::elem (INT i) const
{ 
  if (((unsigned)i) >= length()) error("invalid index");
  return rep->s[i];
}

inline char  String::firstchar() const
{ 
  return elem(0);
}

inline char  String::lastchar() const
{ 
  return elem(length() - 1);
}

// searching

inline INT String::index(char c, INT startpos) const
{
  return search(startpos, length(), c);
}

inline INT String::index(const char* t, INT startpos) const
{   
  return search(startpos, length(), t);
}

inline INT String::index(const String& y, INT startpos) const
{   
  return search(startpos, length(), y.chars(), y.length());
}

inline INT String::index(const SubString& y, INT startpos) const
{   
  return search(startpos, length(), y.chars(), y.length());
}

inline INT String::index(const Regex& r, INT startpos) const
{
  INT unused;  return r.search(chars(), length(), unused, startpos);
}

inline INT String::contains(char c) const
{
  return search(0, length(), c) >= 0;
}

inline INT String::contains(const char* t) const
{   
  return search(0, length(), t) >= 0;
}

inline INT String::contains(const String& y) const
{   
  return search(0, length(), y.chars(), y.length()) >= 0;
}

inline INT String::contains(const SubString& y) const
{   
  return search(0, length(), y.chars(), y.length()) >= 0;
}

inline INT String::contains(char c, INT p) const
{
  return match(p, length(), 0, &c, 1) >= 0;
}

inline INT String::contains(const char* t, INT p) const
{
  return match(p, length(), 0, t) >= 0;
}

inline INT String::contains(const String& y, INT p) const
{
  return match(p, length(), 0, y.chars(), y.length()) >= 0;
}

inline INT String::contains(const SubString& y, INT p) const
{
  return match(p, length(), 0, y.chars(), y.length()) >= 0;
}

inline INT String::contains(const Regex& r) const
{
  INT unused;  return r.search(chars(), length(), unused, 0) >= 0;
}

inline INT String::contains(const Regex& r, INT p) const
{
  return r.match(chars(), length(), p) >= 0;
}


inline INT String::matches(const SubString& y, INT p) const
{
  return match(p, length(), 1, y.chars(), y.length()) >= 0;
}

inline INT String::matches(const String& y, INT p) const
{
  return match(p, length(), 1, y.chars(), y.length()) >= 0;
}

inline INT String::matches(const char* t, INT p) const
{
  return match(p, length(), 1, t) >= 0;
}

inline INT String::matches(char c, INT p) const
{
  return match(p, length(), 1, &c, 1) >= 0;
}

inline INT String::matches(const Regex& r, INT p) const
{
  INT l = (p < 0)? -p : length() - p;
  return r.match(chars(), length(), p) == l;
}


inline INT SubString::contains(const char* t) const
{   
  return S.search(pos, pos+len, t) >= 0;
}

inline INT SubString::contains(const String& y) const
{   
  return S.search(pos, pos+len, y.chars(), y.length()) >= 0;
}

inline INT SubString::contains(const SubString&  y) const
{   
  return S.search(pos, pos+len, y.chars(), y.length()) >= 0;
}

inline INT SubString::contains(char c) const
{
  return S.search(pos, pos+len, c) >= 0;
}

inline INT SubString::contains(const Regex& r) const
{
  INT unused;  return r.search(chars(), len, unused, 0) >= 0;
}

inline INT SubString::matches(const Regex& r) const
{
  return r.match(chars(), len, 0) == len;
}


inline INT String::gsub(const String& pat, const String& r)
{
  return _gsub(pat.chars(), pat.length(), r.chars(), r.length());
}

inline INT String::gsub(const SubString&  pat, const String& r)
{
  return _gsub(pat.chars(), pat.length(), r.chars(), r.length());
}

inline INT String::gsub(const Regex& pat, const String& r)
{
  return _gsub(pat, r.chars(), r.length());
}

inline INT String::gsub(const char* pat, const String& r)
{
  return _gsub(pat, -1, r.chars(), r.length());
}

inline INT String::gsub(const char* pat, const char* r)
{
  return _gsub(pat, -1, r, -1);
}


inline  std::ostream& operator<<(std::ostream& s, const String& x)
{
   s << x.chars(); return s;
}

// a zillion comparison operators

inline INT operator==(const String& x, const String& y) 
{
  return compare(x, y) == 0; 
}

inline INT operator!=(const String& x, const String& y)
{
  return compare(x, y) != 0; 
}

inline INT operator>(const String& x, const String& y)
{
  return compare(x, y) > 0; 
}

inline INT operator>=(const String& x, const String& y)
{
  return compare(x, y) >= 0; 
}

inline INT operator<(const String& x, const String& y)
{
  return compare(x, y) < 0; 
}

inline INT operator<=(const String& x, const String& y)
{
  return compare(x, y) <= 0; 
}

inline INT operator==(const String& x, const SubString&  y) 
{
  return compare(x, y) == 0; 
}

inline INT operator!=(const String& x, const SubString&  y)
{
  return compare(x, y) != 0; 
}

inline INT operator>(const String& x, const SubString&  y)      
{
  return compare(x, y) > 0; 
}

inline INT operator>=(const String& x, const SubString&  y)
{
  return compare(x, y) >= 0; 
}

inline INT operator<(const String& x, const SubString&  y) 
{
  return compare(x, y) < 0; 
}

inline INT operator<=(const String& x, const SubString&  y)
{
  return compare(x, y) <= 0; 
}

inline INT operator==(const String& x, const char* t) 
{
  return compare(x, t) == 0; 
}

inline INT operator!=(const String& x, const char* t) 
{
  return compare(x, t) != 0; 
}

inline INT operator>(const String& x, const char* t)  
{
  return compare(x, t) > 0; 
}

inline INT operator>=(const String& x, const char* t) 
{
  return compare(x, t) >= 0; 
}

inline INT operator<(const String& x, const char* t)  
{
  return compare(x, t) < 0; 
}

inline INT operator<=(const String& x, const char* t) 
{
  return compare(x, t) <= 0; 
}

inline INT operator==(const SubString& x, const String& y) 
{
  return compare(y, x) == 0; 
}

inline INT operator!=(const SubString& x, const String& y)
{
  return compare(y, x) != 0;
}

inline INT operator>(const SubString& x, const String& y)      
{
  return compare(y, x) < 0;
}

inline INT operator>=(const SubString& x, const String& y)     
{
  return compare(y, x) <= 0;
}

inline INT operator<(const SubString& x, const String& y)      
{
  return compare(y, x) > 0;
}

inline INT operator<=(const SubString& x, const String& y)     
{
  return compare(y, x) >= 0;
}

inline INT operator==(const SubString& x, const SubString&  y) 
{
  return compare(x, y) == 0; 
}

inline INT operator!=(const SubString& x, const SubString&  y)
{
  return compare(x, y) != 0;
}

inline INT operator>(const SubString& x, const SubString&  y)      
{
  return compare(x, y) > 0;
}

inline INT operator>=(const SubString& x, const SubString&  y)
{
  return compare(x, y) >= 0;
}

inline INT operator<(const SubString& x, const SubString&  y) 
{
  return compare(x, y) < 0;
}

inline INT operator<=(const SubString& x, const SubString&  y)
{
  return compare(x, y) <= 0;
}

inline INT operator==(const SubString& x, const char* t) 
{
  return compare(x, t) == 0; 
}

inline INT operator!=(const SubString& x, const char* t) 
{
  return compare(x, t) != 0;
}

inline INT operator>(const SubString& x, const char* t)  
{
  return compare(x, t) > 0; 
}

inline INT operator>=(const SubString& x, const char* t) 
{
  return compare(x, t) >= 0; 
}

inline INT operator<(const SubString& x, const char* t)  
{
  return compare(x, t) < 0; 
}

inline INT operator<=(const SubString& x, const char* t) 
{
  return compare(x, t) <= 0; 
}


// a helper needed by at, before, etc.

inline SubString String::_substr(INT first, INT l)
{
  if (first < 0 || (unsigned)(first + l) > length() )
    return SubString(_nilString, 0, 0) ;
  else 
    return SubString(*this, first, l);
}

#endif
