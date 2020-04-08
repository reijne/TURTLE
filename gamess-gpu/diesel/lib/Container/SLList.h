// This may look like C code, but it is really -*- C++ -*-
// WARNING: This file is obsolete.  Use ../SLList.h, if you can.
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


#ifndef _SLList_h
#define _SLList_h 1

#include "../../config.h"

#include "Pix.h"
#include "ListType.defs.h"

#ifndef _SLListNode_h
#define _SLListNode_h 1

template <class ListType>
struct SLListNode
{
  SLListNode<ListType>*         tl;
  ListType                    hd;
                         SLListNode() { }
                         SLListNode(const ListType& h, SLListNode<ListType>* t = 0);
                         ~SLListNode() { }
};


template <class ListType>
inline SLListNode<ListType>::SLListNode(const ListType& h, SLListNode<ListType>* t)
:tl(t), hd(h) {}

//typedef SLListNode* SLListNodePtr;

#endif


template <class ListType>
class SLList
{
protected:
  SLListNode<ListType>*        last;

public:
                        SLList();
                        SLList(const SLList<ListType>& a);
                        virtual ~SLList();

  SLList&            operator = (const SLList<ListType>& a);

  INT                   empty() const;
  INT                   length() const;

  void                  clear();

  Pix                   prepend(ListType& item);
  Pix                   append(ListType& item);

  void                  join(SLList&);

  Pix                   prepend(SLListNode<ListType>*);
  Pix                   append(SLListNode<ListType>*);

  ListType&                  operator () (Pix p);
  ListType&                  operator () (Pix p) const;
  Pix                   first() const;
  void                  next(Pix& p) const;
  INT                   owns(Pix p) const;
  Pix                   ins_after(Pix p, ListType& item);
  void                  del_after(Pix p);

  ListType&                  front() const;
  ListType&                  rear() const;
  ListType                   remove_front();
  INT                   remove_front(ListType& x);
  void                  del_front();

  void                  error(const char* msg) const;
  INT                   OK();
};

template <class ListType>
inline SLList<ListType>::~SLList()
{
  clear();
}

template <class ListType>
inline SLList<ListType>::SLList()
{
  last = 0;
}

template <class ListType>
inline INT SLList<ListType>::empty() const
{
  return last == 0;
}


template <class ListType>
inline Pix SLList<ListType>::first() const
{
  return (last == 0)? 0 : Pix(last->tl);
}

template <class ListType>
inline void SLList<ListType>::next(Pix& p) const
{
  p = (p == 0 || p == last)? 0 : Pix(((SLListNode<ListType>*)(p))->tl);
}

template <class ListType>
inline ListType& SLList<ListType>::operator () (Pix p)
{
  if (p == 0) error("null Pix");
  return ((SLListNode<ListType>*)(p))->hd;
}

template <class ListType>
inline ListType& SLList<ListType>::operator () (Pix p) const
{
  if (p == 0) error("null Pix");
  return ((SLListNode<ListType>*)(p))->hd;
}

template <class ListType>
inline ListType& SLList<ListType>::front() const
{
  if (last == 0) error("front: empty list");
  return last->tl->hd;
}

template <class ListType>
inline ListType& SLList<ListType>::rear() const
{
  if (last == 0) error("rear: empty list");
  return last->hd;
}

#endif
