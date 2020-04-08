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

#include "builtin.h"
#include "Map.h"


template <class KeyType, class ObjectType>
Pix Map<KeyType, ObjectType>::seek(KeyType  item)
{
  Pix i;
  for (i = first(); i != 0 && !(KeyTypeEQ(key(i), item)); next(i));
  return i;
}

template <class KeyType, class ObjectType>
INT Map<KeyType, ObjectType>::owns(Pix idx)
{
  if (idx == 0) return 0;
  for (Pix i = first(); i; next(i)) if (i == idx) return 1;
  return 0;
}

template <class KeyType, class ObjectType>
void Map<KeyType, ObjectType>::clear()
{
  Pix i = first(); 
  while (i != 0)
  {
    del(key(i));
    i = first();
  }
}

template <class KeyType, class ObjectType>
INT Map<KeyType, ObjectType>::contains (KeyType  item)
{
  return seek(item) != 0;
}


template <class KeyType, class ObjectType>
void Map<KeyType, ObjectType>::error(const char* msg)
{
  (*lib_error_handler)("Map", msg);
}

/*
#ifdef __GNUC__

#include "TableKey.h"
#include "TAccessList.h"

template class Map<TableKey, TAccessList>;



#include "ConfigurationMOType.h"
class EnergyEntry;
template class Map<ConfigurationMOType, EnergyEntry *>;



#endif
*/

