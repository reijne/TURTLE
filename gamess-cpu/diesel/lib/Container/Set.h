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
Foundation, 675 Mass Ave, Cambridge, MA 02139, USA.
*/


#ifndef _Set_h
#define _Set_h 1

#include "../../config.h"

#include "Pix.h"
#include "SetType.defs.h"


template <class SetType>
class Set
{
protected:

  INT                   count;

public:
  virtual              ~Set();


  INT                   length() const;                // current number of items
  INT                   empty() const;

  virtual Pix           add(SetType& item) = 0;      // add item; return Pix
  virtual void          del(SetType& item) = 0;      // delete item
  virtual INT           contains(SetType& item);     // is item in set?

  virtual void          clear();                 // delete all items

  virtual Pix           first() const = 0;             // Pix of first item or 0
  virtual void          next(Pix& i) const = 0;        // advance to next or 0
  virtual SetType&          operator () (Pix i) const = 0; // access item at i

  virtual INT           owns(Pix i);             // is i a valid Pix  ?
  virtual Pix           seek(SetType& item) const;         // Pix of item

  void                  operator |= (Set& b); // add all items in b
  void                  operator -= (Set& b); // delete items also in b
  void                  operator &= (Set& b); // delete items not in b

  INT                   operator == (Set& b);
  INT                   operator != (Set& b);
  INT                   operator <= (Set& b); 

  void                  error(const char* msg) const;
  virtual INT           OK() = 0;                // rep invariant
};

template <class SetType>
inline Set<SetType>::~Set() {}

template <class SetType>
inline INT Set<SetType>::length() const
{
  return count;
}

template <class SetType>
inline INT Set<SetType>::empty() const
{
  return count == 0;
}

#endif
