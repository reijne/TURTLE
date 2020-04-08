// This may look like C code, but it is really -*- C++ -*-
// WARNING: This file is obsolete.  Use ../SLList.cc, if you can.
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

#include <limits.h>
#include <iostream>

#include "builtin.h"
#include "SLList.h"

using namespace std;

template <class ListType>
void SLList<ListType>::error(const char* msg) const
{
  (*lib_error_handler)("SLList", msg);
}

template <class ListType>
INT SLList<ListType>::length() const
{
  INT l = 0;
  SLListNode<ListType>* t = last;
  if (t != 0) do { ++l; t = t->tl; } while (t != last);
  return l;
}

template <class ListType>
SLList<ListType>::SLList(const SLList<ListType>& a)
{
  if (a.last == 0)
    last = 0;
  else
  {
    SLListNode<ListType>* p = a.last->tl;
    SLListNode<ListType>* h = new SLListNode<ListType>(p->hd);
    last = h;
    for (;;)
    {
      if (p == a.last)
      {
        last->tl = h;
        return;
      }
      p = p->tl;
      SLListNode<ListType>* n = new SLListNode<ListType>(p->hd);
      last->tl = n;
      last = n;
    }
  }
}

template <class ListType>
SLList<ListType>& SLList<ListType>::operator = (const SLList<ListType>& a)
{
  if (last != a.last)
  {
    clear();
    if (a.last != 0)
    {
      SLListNode<ListType>* p = a.last->tl;
      SLListNode<ListType>* h = new SLListNode<ListType>(p->hd);
      last = h;
      for (;;)
      {
        if (p == a.last)
        {
          last->tl = h;
          break;
        }
        p = p->tl;
        SLListNode<ListType>* n = new SLListNode<ListType>(p->hd);
        last->tl = n;
        last = n;
      }
    }
  }
  return *this;
}

template <class ListType>
void SLList<ListType>::clear()
{
  if (last == 0)
    return;

  SLListNode<ListType>* p = last->tl;
  last->tl = 0;
  last = 0;

  while (p != 0)
  {
    SLListNode<ListType>* nxt = p->tl;
    delete(p);
    p = nxt;
  }
}


template <class ListType>
Pix SLList<ListType>::prepend(ListType& item)
{
  SLListNode<ListType>* t = new SLListNode<ListType>(item);
  if (last == 0)
    t->tl = last = t;
  else
  {
    t->tl = last->tl;
    last->tl = t;
  }
  return Pix(t);
}


template <class ListType>
Pix SLList<ListType>::prepend(SLListNode<ListType>* t)
{
  if (t == 0) return 0;
  if (last == 0)
    t->tl = last = t;
  else
  {
    t->tl = last->tl;
    last->tl = t;
  }
  return Pix(t);
}


template <class ListType>
Pix SLList<ListType>::append(ListType& item)
{
  SLListNode<ListType>* t = new SLListNode<ListType>(item);
  if (last == 0)
    t->tl = last = t;
  else
  {
    t->tl = last->tl;
    last->tl = t;
    last = t;
  }
  return Pix(t);
}

template <class ListType>
Pix SLList<ListType>::append(SLListNode<ListType>* t)
{
  if (t == 0) return 0;
  if (last == 0)
    t->tl = last = t;
  else
  {
    t->tl = last->tl;
    last->tl = t;
    last = t;
  }
  return Pix(t);
}

template <class ListType>
void SLList<ListType>::join(SLList<ListType>& b)
{
  SLListNode<ListType>* t = b.last;
  b.last = 0;
  if (last == 0)
    last = t;
  else if (t != 0)
  {
    SLListNode<ListType>* f = last->tl;
    last->tl = t->tl;
    t->tl = f;
    last = t;
  }
}

template <class ListType>
Pix SLList<ListType>::ins_after(Pix p, ListType& item)
{
  SLListNode<ListType>* u = (SLListNode<ListType>*)p;
  SLListNode<ListType>* t = new SLListNode<ListType>(item);
  if (last == 0)
    t->tl = last = t;
  else if (u == 0) // ins_after 0 means prepend
  {
    t->tl = last->tl;
    last->tl = t;
  }
  else
  {
    t->tl = u->tl;
    u->tl = t;
    if (u == last) 
      last = t;
  }
  return Pix(t);
}


template <class ListType>
void SLList<ListType>::del_after(Pix p)
{
  SLListNode<ListType>* u = (SLListNode<ListType>*)p;
  if (last == 0 || u == last) error("cannot del_after last");
  if (u == 0) u = last; // del_after 0 means delete first
  SLListNode<ListType>* t = u->tl;
  if (u == t)
    last = 0;
  else
  {
    u->tl = t->tl;
    if (last == t)
      last = u;
  }
  delete t;
}

template <class ListType>
INT SLList<ListType>::owns(Pix p) const
{
  SLListNode<ListType>* t = last;
  if (t != 0 && p != 0)
  {
    do
    {
      if (Pix(t) == p) return 1;
      t = t->tl;
    } while (t != last);
  }
  return 0;
}

template <class ListType>
ListType SLList<ListType>::remove_front()
{
  if (last == 0) error("remove_front of empty list");
  SLListNode<ListType>* t = last->tl;
  ListType res = t->hd;
  if (t == last)
    last = 0;
  else
    last->tl = t->tl;
  delete t;
  return res;
}

template <class ListType>
INT SLList<ListType>::remove_front(ListType& x)
{
  if (last == 0)
    return 0;
  else
  {
    SLListNode<ListType>* t = last->tl;
    x = t->hd;
    if (t == last)
      last = 0;
    else
      last->tl = t->tl;
    delete t;
    return 1;
  }
}


template <class ListType>
void SLList<ListType>::del_front()
{
  if (last == 0) error("del_front of empty list");
  SLListNode<ListType>* t = last->tl;
  if (t == last)
    last = 0;
  else
    last->tl = t->tl;
  delete t;
}

template <class ListType>
INT SLList<ListType>::OK()
{
  INT v = 1;
  if (last != 0)
  {
    SLListNode<ListType>* t = last;
    long count = LONG_MAX;      // Lots of chances to find last!
    do
    {
      count--;
      t = t->tl;
    } while (count > 0 && t != last);
    v &= count > 0;
  }
  if (!v) error("invariant failure");
  return v;
}





/*
#include "extEntry.h"
#include "MatchingMOListIterator.h"

template class SLList<extEntry *>;
template class SLList<MatchingMOList>;


#include "../../../GUI/vis/calc/SpaceGridSurfacePointList.h"
#include "../../../GUI/vis/calc/SurfacePoint.h"
template class SLList<SpaceGridSurfacePointList *>;
template class SLList<SurfacePoint *>;
*/
