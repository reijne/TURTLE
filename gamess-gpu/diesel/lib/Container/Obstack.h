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


#ifndef _Obstack_h
#define _Obstack_h 1
#include <string>
#include "../../config.h"

#undef OK

class Obstack
{
  struct _obstack_chunk
  {
    char*           limit;
    _obstack_chunk* prev;
    char            contents[4];
  };

protected:
  LONG_INT	          chunksize;
  _obstack_chunk* chunk;
  char*	          objectbase;
  char*	          nextfree;
  char*	          chunklimit;
  INT             alignmentmask;

  void  _free(void* obj);
  void  newchunk(INT size);

public:
        Obstack(INT size = 4080, INT alignment = 4); // 4080=4096-mallocslop

        ~Obstack();

  void* base();
  void* next_free();
  INT   alignment_mask();
  INT   chunk_size();
  INT   size();
  INT   room();
  INT   contains(void* p);      // does Obstack hold pointer p?

  void  grow(const void* data, INT size);
  void  grow(const void* data, INT size, char terminator);
  void  grow(const char* s);
  void  grow(char c);
  void  grow_fast(char c);
  void  blank(INT size);
  void  blank_fast(INT size);

  void* finish();
  void* finish(char terminator);

  void* copy(const void* data, INT size);
  void* copy(const void* data, INT size, char terminator);
  void* copy(const char* s);
  void* copy(char c);
  void* alloc(INT size);

  void  free(void* obj);
  void  shrink(INT size = 1); // suggested by ken@cs.rochester.edu

  INT   OK();                   // rep invariant
};


inline Obstack::~Obstack()
{
  _free(0); 
}

inline void* Obstack::base()
{
  return objectbase; 
}

inline void* Obstack::next_free()
{
  return nextfree; 
}

inline INT Obstack::alignment_mask()
{
  return alignmentmask; 
}

inline INT Obstack::chunk_size()
{
  return chunksize; 
}

inline INT Obstack::size()
{
  return nextfree - objectbase; 
}

inline INT Obstack::room()
{
  return chunklimit - nextfree; 
}

inline void Obstack:: grow(const void* data, INT size)
{
  if (nextfree+size > chunklimit) 
    newchunk(size);
  memcpy(nextfree, data, size);
  nextfree += size; 
}

inline void Obstack:: grow(const void* data, INT size, char terminator)
{
  if (nextfree+size+1 > chunklimit) 
    newchunk(size+1);
  memcpy(nextfree, data, size);
  nextfree += size; 
  *(nextfree)++ = terminator; 
}

inline void Obstack:: grow(const char* s)
{
  grow((const void*)s, strlen(s), 0); 
}

inline void Obstack:: grow(char c)
{
  if (nextfree+1 > chunklimit) 
    newchunk(1); 
  *(nextfree)++ = c; 
}

inline void Obstack:: blank(INT size)
{
  if (nextfree+size > chunklimit) 
    newchunk(size);
  nextfree += size; 
}

inline void* Obstack::finish(char terminator)
{
  grow(terminator); 
  return finish(); 
}

inline void* Obstack::copy(const void* data, INT size)
{
  grow (data, size);
  return finish(); 
}

inline void* Obstack::copy(const void* data, INT size, char terminator)
{
  grow(data, size, terminator); 
  return finish(); 
}

inline void* Obstack::copy(const char* s)
{
  grow((const void*)s, strlen(s), 0); 
  return finish(); 
}

inline void* Obstack::copy(char c)
{
  grow(c);
  return finish(); 
}

inline void* Obstack::alloc(INT size)
{
  blank(size);
  return finish(); 
}

inline void Obstack:: free(void* obj)     
{
  if (obj >= (void*)chunk && obj<(void*)chunklimit)
    nextfree = objectbase = (char *) obj;
  else 
    _free(obj); 
}

inline void Obstack:: grow_fast(char c)
{
  *(nextfree)++ = c; 
}

inline void Obstack:: blank_fast(INT size)
{
  nextfree += size; 
}

inline void Obstack:: shrink(INT size) // from ken@cs.rochester.edu
{
  if (nextfree >= objectbase + size)
    nextfree -= size;
}

#endif
