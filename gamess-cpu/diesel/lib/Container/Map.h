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


#ifndef _Map_h
#define _Map_h 1

#include "../../config.h"

#include "Pix.h"
#include "KeyType.defs.h"

template <class KeyType, class ObjectType>
class Map
{
protected:
  INT                   count;
  ObjectType            def;

public:
                        Map(ObjectType  dflt);
  virtual              ~Map();

  INT                   length() const; // current number of items
  INT                   empty() const;

  virtual INT           contains(KeyType  key);      // is key mapped?

  virtual void          clear();                 // delete all items

  virtual ObjectType&   operator [] (KeyType  key) = 0; // access contents by key

  virtual void          del(KeyType  key) = 0;       // delete entry

  virtual Pix           first() = 0;             // Pix of first item or 0
  virtual void          next(Pix& i) = 0;        // advance to next or 0
  virtual KeyType&      key(Pix i) = 0;          // access key at i
  virtual ObjectType&   contents(Pix i) = 0;     // access contents at i

  virtual INT           owns(Pix i);             // is i a valid Pix  ?
  virtual Pix           seek(KeyType  key);          // Pix of key

  ObjectType&           dflt();                  // access default val

  void                  error(const char* msg);
  virtual INT           OK() = 0;                // rep invariant
};


template <class KeyType, class ObjectType>
inline Map<KeyType, ObjectType>::~Map() {}

template <class KeyType, class ObjectType>
inline INT Map<KeyType, ObjectType>::length() const
{
  return count;
}

template <class KeyType, class ObjectType>
inline INT Map<KeyType, ObjectType>::empty() const
{
  return count == 0;
}

template <class KeyType, class ObjectType>
inline ObjectType& Map<KeyType, ObjectType>::dflt()
{
  return def;
}

template <class KeyType, class ObjectType>
inline Map<KeyType, ObjectType>::Map(ObjectType  dflt) :def(dflt)
{
  count = 0;
}


#endif
