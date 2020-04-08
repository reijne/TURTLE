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


#ifndef _SetTypeAVL_h
#define _SetTypeAVL_h 1

#include "../../config.h"

#include "Set.h"

template <class SetType>
struct AVLSetNode
{
  AVLSetNode<SetType>*         lt;
  AVLSetNode<SetType>*         rt;
  SetType                 item;
  char                stat;
                      AVLSetNode(SetType& h, AVLSetNode<SetType>* l=0, AVLSetNode<SetType>* r=0);
                      AVLSetNode(SetType& h, INT dummy);
                      ~AVLSetNode();
};



template <class SetType>
inline AVLSetNode<SetType>::AVLSetNode(SetType& h, AVLSetNode<SetType>* l, AVLSetNode<SetType>* r)
:lt(l), rt(r), item(h), stat(0) {}


template <class SetType>
inline AVLSetNode<SetType>::~AVLSetNode() {}

//typedef AVLSetNode<SetType>* AVLSetNodePtr;


template <class SetType>
class AVLSet : public Set<SetType>
{
private:
SetType*   _target_item;     // add/del_item target
AVLSetNode<SetType>* _found_node; // returned added/deleted node
AVLSetNode<SetType>** _hold_nodes;       // used for rebuilding trees


AVLSetNode<SetType>* _treeify(INT n);
AVLSetNode<SetType>* _do_treeify(INT lo, INT hi, INT& h);

protected:
  AVLSetNode<SetType>*   root;

                AVLSet(AVLSetNode<SetType>* p, INT l);

  AVLSetNode<SetType>*   leftmost() const;
  AVLSetNode<SetType>*   rightmost() const;
  AVLSetNode<SetType>*   pred(AVLSetNode<SetType>* t) const;
  AVLSetNode<SetType>*   succ(AVLSetNode<SetType>* t) const;
  void          _kill(AVLSetNode<SetType>* t);
  void          _add(AVLSetNode<SetType>*& t);
  void          _del(AVLSetNode<SetType>* p, AVLSetNode<SetType>*& t);

public:
                AVLSet();
                AVLSet(AVLSet<SetType>& a);
                ~AVLSet();
				// "virtual" constructor:
				virtual AVLSet<SetType> *	new_AVLSet() { return new AVLSet<SetType>(); }

  Pix           add(SetType& item);
  void          del(SetType& item);
  INT           contains(SetType& item);

  void          clear();

  Pix           first() const;
  void          next(Pix& i) const;
  SetType&          operator () (Pix i) const;
  INT           owns(Pix i);
  Pix           seek(SetType & item) const;

  Pix           last() const;
  void          prev(Pix& i);

  void          operator |= (AVLSet& b);
  void          operator -= (AVLSet& b);
  void          operator &= (AVLSet& b);

  INT           operator == (AVLSet& b);
  INT           operator != (AVLSet& b);
  INT           operator <= (AVLSet& b); 

  INT           OK();
};

template <class SetType>
inline AVLSet<SetType>::~AVLSet()
{
  _kill(root);
}

template <class SetType>
inline AVLSet<SetType>::AVLSet()
{
  root = 0;
  this->count = 0;
}

template <class SetType>
inline AVLSet<SetType>::AVLSet(AVLSetNode<SetType>* p, INT l)
{
  root = p;
  this->count = l;
}

template <class SetType>
inline INT AVLSet<SetType>::operator != (AVLSet& b)
{
  return ! ((*this) == b);
}

template <class SetType>
inline Pix AVLSet<SetType>::first() const
{
  return Pix(leftmost());
}

template <class SetType>
inline Pix AVLSet<SetType>::last() const
{
  return Pix(rightmost());
}

template <class SetType>
inline void AVLSet<SetType>::next(Pix& i) const
{
  if (i != 0) i = Pix(succ((AVLSetNode<SetType>*)i));
}

template <class SetType>
inline void AVLSet<SetType>::prev(Pix& i)
{
  if (i != 0) i = Pix(pred((AVLSetNode<SetType>*)i));
}

template <class SetType>
inline SetType& AVLSet<SetType>::operator () (Pix i) const
{
  if (i == 0) this->error("null Pix");
  return ((AVLSetNode<SetType>*)i)->item;
}

template <class SetType>
inline void AVLSet<SetType>::clear()
{
  _kill(root);
  this->count = 0;
  root = 0;
}

template <class SetType>
inline INT AVLSet<SetType>::contains(SetType& key)
{
  return seek(key) != 0;
}

#endif
