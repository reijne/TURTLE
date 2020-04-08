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

#include <iostream>

using namespace std;

#include "AVLSet.h"
#include <stdlib.h>


/*
 constants & inlines for maintaining balance & thread status in tree nodes
*/

#define AVLBALANCEMASK    3
#define AVLBALANCED       0
#define AVLLEFTHEAVY      1
#define AVLRIGHTHEAVY     2

#define LTHREADBIT        4
#define RTHREADBIT        8


template <class SetType>
static inline INT bf(AVLSetNode<SetType>* t)
{
  return t->stat & AVLBALANCEMASK;
}

template <class SetType>
static inline void set_bf(AVLSetNode<SetType>* t, INT b)
{
  t->stat = (t->stat & ~AVLBALANCEMASK) | (b & AVLBALANCEMASK);
}


template <class SetType>
static inline INT rthread(AVLSetNode<SetType>* t)
{
  return t->stat & RTHREADBIT;
}

template <class SetType>
static inline void set_rthread(AVLSetNode<SetType>* t, INT b)
{
  if (b)
    t->stat |= RTHREADBIT;
  else
    t->stat &= ~RTHREADBIT;
}

template <class SetType>
static inline INT lthread(AVLSetNode<SetType>* t)
{
  return t->stat & LTHREADBIT;
}

template <class SetType>
static inline void set_lthread(AVLSetNode<SetType>* t, INT b)
{
  if (b)
    t->stat |= LTHREADBIT;
  else
    t->stat &= ~LTHREADBIT;
}

/*
 traversal primitives
*/


template <class SetType>
AVLSetNode<SetType>* AVLSet<SetType>::leftmost() const
{
  AVLSetNode<SetType>* t = root;
  if (t != 0) while (t->lt != 0) t = t->lt;
  return t;
}

template <class SetType>
AVLSetNode<SetType>* AVLSet<SetType>::rightmost() const
{
  AVLSetNode<SetType>* t = root;
  if (t != 0) while (t->rt != 0) t = t->rt;
  return t;
}

template <class SetType>
AVLSetNode<SetType>* AVLSet<SetType>::succ(AVLSetNode<SetType>* t) const
{
  AVLSetNode<SetType>* r = t->rt;
  if (!rthread(t)) while (!lthread(r)) r = r->lt;
  return r;
}

template <class SetType>
AVLSetNode<SetType>* AVLSet<SetType>::pred(AVLSetNode<SetType>* t) const
{
  AVLSetNode<SetType>* l = t->lt;
  if (!lthread(t)) while (!rthread(l)) l = l->rt;
  return l;
}


template <class SetType>
Pix AVLSet<SetType>::seek(SetType& key) const
{
  AVLSetNode<SetType>* t = root;
  if (t == 0)
    return 0;
  for (;;)
  {
    INT cmp = SetTypeCMP(key, t->item);
    if (cmp == 0)
      return Pix(t);
    else if (cmp < 0)
    {
      if (lthread(t))
        return 0;
      else
        t = t->lt;
    }
    else if (rthread(t))
      return 0;
    else
      t = t->rt;
  }
}


/*
 The combination of threads and AVL bits make adding & deleting
 interesting, but very awkward.

 We use the following statics to avoid passing them around recursively
*/

static INT _need_rebalancing;   // to send back balance info from rec. calls
static INT    _already_found;   // for deletion subcases

static INT  _max_hold_index;              // # elements-1 in _hold_nodes


template <class SetType>
void AVLSet<SetType>:: _add(AVLSetNode<SetType>*& t)
{
  INT cmp = SetTypeCMP(*_target_item, t->item);
  if (cmp == 0)
  {
    _found_node = t;
    return;
  }
  else if (cmp < 0)
  {
    if (lthread(t))
    {
      ++this->count;
      _found_node = new AVLSetNode<SetType>(*_target_item);
      set_lthread(_found_node, 1);
      set_rthread(_found_node, 1);
      _found_node->lt = t->lt;
      _found_node->rt = t;
      t->lt = _found_node;
      set_lthread(t, 0);
      _need_rebalancing = 1;
    }
    else
      _add(t->lt);
    if (_need_rebalancing)
    {
      switch(bf(t))
      {
      case AVLRIGHTHEAVY:
        set_bf(t, AVLBALANCED);
        _need_rebalancing = 0;
        return;
      case AVLBALANCED:
        set_bf(t, AVLLEFTHEAVY);
        return;
      case AVLLEFTHEAVY:
	{
        AVLSetNode<SetType>* l = t->lt;
        if (bf(l) == AVLLEFTHEAVY)
        {
          if (rthread(l))
            t->lt = l;
          else
            t->lt = l->rt;
          set_lthread(t, rthread(l));
          l->rt = t;
          set_rthread(l, 0);
          set_bf(t, AVLBALANCED);
          set_bf(l, AVLBALANCED);
          t = l;
          _need_rebalancing = 0;
        }
        else
        {
          AVLSetNode<SetType>* r = l->rt;
          set_rthread(l, lthread(r));
          if (lthread(r))
            l->rt = r;
          else
            l->rt = r->lt;
          r->lt = l;
          set_lthread(r, 0);
          set_lthread(t, rthread(r));
          if (rthread(r))
            t->lt = r;
          else
            t->lt = r->rt;
          r->rt = t;
          set_rthread(r, 0);
          if (bf(r) == AVLLEFTHEAVY)
            set_bf(t, AVLRIGHTHEAVY);
          else
            set_bf(t, AVLBALANCED);
          if (bf(r) == AVLRIGHTHEAVY)
            set_bf(l, AVLLEFTHEAVY);
          else
            set_bf(l, AVLBALANCED);
          set_bf(r, AVLBALANCED);
          t = r;
          _need_rebalancing = 0;
          return;
        }
	}
      }
    }
  }
  else
  {
    if (rthread(t))
    {
      ++this->count;
      _found_node = new AVLSetNode<SetType>(*_target_item);
      set_rthread(t, 0);
      set_lthread(_found_node, 1);
      set_rthread(_found_node, 1);
      _found_node->lt = t;
      _found_node->rt = t->rt;
      t->rt = _found_node;
      _need_rebalancing = 1;
    }
    else
      _add(t->rt);
    if (_need_rebalancing)
    {
      switch(bf(t))
      {
      case AVLLEFTHEAVY:
        set_bf(t, AVLBALANCED);
        _need_rebalancing = 0;
        return;
      case AVLBALANCED:
        set_bf(t, AVLRIGHTHEAVY);
        return;
      case AVLRIGHTHEAVY:
	{
        AVLSetNode<SetType>* r = t->rt;
        if (bf(r) == AVLRIGHTHEAVY)
        {
          if (lthread(r))
            t->rt = r;
          else
            t->rt = r->lt;
          set_rthread(t, lthread(r));
          r->lt = t;
          set_lthread(r, 0);
          set_bf(t, AVLBALANCED);
          set_bf(r, AVLBALANCED);
          t = r;
          _need_rebalancing = 0;
        }
        else
        {
          AVLSetNode<SetType>* l = r->lt;
          set_lthread(r, rthread(l));
          if (rthread(l))
            r->lt = l;
          else
            r->lt = l->rt;
          l->rt = r;
          set_rthread(l, 0);
          set_rthread(t, lthread(l));
          if (lthread(l))
            t->rt = l;
          else
            t->rt = l->lt;
          l->lt = t;
          set_lthread(l, 0);
          if (bf(l) == AVLRIGHTHEAVY)
            set_bf(t, AVLLEFTHEAVY);
          else
            set_bf(t, AVLBALANCED);
          if (bf(l) == AVLLEFTHEAVY)
            set_bf(r, AVLRIGHTHEAVY);
          else
            set_bf(r, AVLBALANCED);
          set_bf(l, AVLBALANCED);
          t = l;
          _need_rebalancing = 0;
          return;
        }
	}
      }
    }
  }
}

    
template <class SetType>
Pix AVLSet<SetType>::add(SetType& item)
{
  if (root == 0)
  {
    ++this->count;
    root = new AVLSetNode<SetType>(item);
    set_rthread(root, 1);
    set_lthread(root, 1);
    return Pix(root);
  }
  else
  {
    _target_item = &item;
    _need_rebalancing = 0;
    _add(root);
    return Pix(_found_node);
  }
}


template <class SetType>
void AVLSet<SetType>::_del(AVLSetNode<SetType>* par, AVLSetNode<SetType>*& t)
{
  INT comp;
  if (_already_found)
  {
    if (rthread(t))
      comp = 0;
    else
      comp = 1;
  }
  else 
    comp = SetTypeCMP(*_target_item, t->item);
  if (comp == 0)
  {
    if (lthread(t) && rthread(t))
    {
      _found_node = t;
      if (t == par->lt)
      {
        set_lthread(par, 1);
        par->lt = t->lt;
      }
      else
      {
        set_rthread(par, 1);
        par->rt = t->rt;
      }
      _need_rebalancing = 1;
      return;
    }
    else if (lthread(t))
    {
      _found_node = t;
      AVLSetNode<SetType>* s = succ(t);
      if (s != 0 && lthread(s))
        s->lt = t->lt;
      t = t->rt;
      _need_rebalancing = 1;
      return;
    }
    else if (rthread(t))
    {
      _found_node = t;
      AVLSetNode<SetType>* p = pred(t);
      if (p != 0 && rthread(p))
        p->rt = t->rt;
      t = t->lt;
      _need_rebalancing = 1;
      return;
    }
    else                        // replace item & find someone deletable
    {
      AVLSetNode<SetType>* p = pred(t);
      t->item = p->item;
      _already_found = 1;
      comp = -1;                // fall through below to left
    }
  }

  if (comp < 0)
  {
    if (lthread(t))
      return;
    _del(t, t->lt);
    if (!_need_rebalancing)
      return;
    switch (bf(t))
    {
    case AVLLEFTHEAVY:
      set_bf(t, AVLBALANCED);
      return;
    case AVLBALANCED:
      set_bf(t, AVLRIGHTHEAVY);
      _need_rebalancing = 0;
      return;
    case AVLRIGHTHEAVY:
      {
      AVLSetNode<SetType>* r = t->rt;
      switch (bf(r))
      {
      case AVLBALANCED:
        if (lthread(r))
          t->rt = r;
        else
          t->rt = r->lt;
        set_rthread(t, lthread(r));
        r->lt = t;
        set_lthread(r, 0);
        set_bf(t, AVLRIGHTHEAVY);
        set_bf(r, AVLLEFTHEAVY);
        _need_rebalancing = 0;
        t = r;
        return;
      case AVLRIGHTHEAVY:
        if (lthread(r))
          t->rt = r;
        else
          t->rt = r->lt;
        set_rthread(t, lthread(r));
        r->lt = t;
        set_lthread(r, 0);
        set_bf(t, AVLBALANCED);
        set_bf(r, AVLBALANCED);
        t = r;
        return;
      case AVLLEFTHEAVY:
	{
        AVLSetNode<SetType>* l = r->lt;
        set_lthread(r, rthread(l));
        if (rthread(l))
          r->lt = l;
        else
          r->lt = l->rt;
        l->rt = r;
        set_rthread(l, 0);
        set_rthread(t, lthread(l));
        if (lthread(l))
          t->rt = l;
        else
          t->rt = l->lt;
        l->lt = t;
        set_lthread(l, 0);
        if (bf(l) == AVLRIGHTHEAVY)
          set_bf(t, AVLLEFTHEAVY);
        else
          set_bf(t, AVLBALANCED);
        if (bf(l) == AVLLEFTHEAVY)
          set_bf(r, AVLRIGHTHEAVY);
        else
          set_bf(r, AVLBALANCED);
        set_bf(l, AVLBALANCED);
        t = l;
        return;
	}
      }
    }
    }
  }
  else
  {
    if (rthread(t))
      return;
    _del(t, t->rt);
    if (!_need_rebalancing)
      return;
    switch (bf(t))
    {
    case AVLRIGHTHEAVY:
      set_bf(t, AVLBALANCED);
      return;
    case AVLBALANCED:
      set_bf(t, AVLLEFTHEAVY);
      _need_rebalancing = 0;
      return;
    case AVLLEFTHEAVY:
      {
      AVLSetNode<SetType>* l = t->lt;
      switch (bf(l))
      {
      case AVLBALANCED:
        if (rthread(l))
          t->lt = l;
        else
          t->lt = l->rt;
        set_lthread(t, rthread(l));
        l->rt = t;
        set_rthread(l, 0);
        set_bf(t, AVLLEFTHEAVY);
        set_bf(l, AVLRIGHTHEAVY);
        _need_rebalancing = 0;
        t = l;
        return;
      case AVLLEFTHEAVY:
        if (rthread(l))
          t->lt = l;
        else
          t->lt = l->rt;
        set_lthread(t, rthread(l));
        l->rt = t;
        set_rthread(l, 0);
        set_bf(t, AVLBALANCED);
        set_bf(l, AVLBALANCED);
        t = l;
        return;
      case AVLRIGHTHEAVY:
	{
        AVLSetNode<SetType>* r = l->rt;
        set_rthread(l, lthread(r));
        if (lthread(r))
          l->rt = r;
        else
          l->rt = r->lt;
        r->lt = l;
        set_lthread(r, 0);
        set_lthread(t, rthread(r));
        if (rthread(r))
          t->lt = r;
        else
          t->lt = r->rt;
        r->rt = t;
        set_rthread(r, 0);
        if (bf(r) == AVLLEFTHEAVY)
          set_bf(t, AVLRIGHTHEAVY);
        else
          set_bf(t, AVLBALANCED);
        if (bf(r) == AVLRIGHTHEAVY)
          set_bf(l, AVLLEFTHEAVY);
        else
          set_bf(l, AVLBALANCED);
        set_bf(r, AVLBALANCED);
        t = r;
        return;
	}
      }
      }
    }
  }
}

        

template <class SetType>
void AVLSet<SetType>::del(SetType& item)
{
  if (root == 0) return;
  _need_rebalancing = 0;
  _already_found = 0;
  _found_node = 0;
  _target_item = &item;
  _del(root, root);
  if (_found_node)
  {
    delete(_found_node);
    if (--this->count == 0)
      root = 0;
  }
}

// build an ordered array of pointers to tree nodes back into a tree
// we know that at least one element exists

template <class SetType>
AVLSetNode<SetType>* AVLSet<SetType>::_do_treeify(INT lo, INT hi, INT& h)
{
  INT lh, rh;
  INT mid = (lo + hi) / 2;
  AVLSetNode<SetType>* t = _hold_nodes[mid];
  if (lo > mid - 1)
  {
    set_lthread(t, 1);
    if (mid == 0)
      t->lt = 0;
    else
      t->lt = _hold_nodes[mid-1];
    lh = 0;
  }
  else
  {
    set_lthread(t, 0);
    t->lt = _do_treeify(lo, mid-1, lh);
  }
  if (hi < mid + 1)
  {
    set_rthread(t, 1);
    if (mid == _max_hold_index)
      t->rt = 0;
    else
      t->rt = _hold_nodes[mid+1];
    rh = 0;
  }
  else 
  {
    set_rthread(t, 0);
    t->rt = _do_treeify(mid+1, hi, rh);
  }
  if (lh == rh)
  {
    set_bf(t, AVLBALANCED);
    h = lh + 1;
  }
  else if (lh == rh - 1)
  {
    set_bf(t, AVLRIGHTHEAVY);
    h = rh + 1;
  }
  else if (rh == lh - 1)
  {
    set_bf(t, AVLLEFTHEAVY);
    h = lh + 1;
  }
  else                          // can't happen
    abort();

  return t;
}

template <class SetType>
AVLSetNode<SetType>* AVLSet<SetType>::_treeify(INT n)
{
  AVLSetNode<SetType>* t;
  if (n == 0)
    t = 0;
  else
  {
    INT b;
    _max_hold_index = n-1;
    t = _do_treeify(0, _max_hold_index, b);
  }
  delete[] _hold_nodes;
  return t;
}


template <class SetType>
void AVLSet<SetType>::_kill(AVLSetNode<SetType>* t)
{
  if (t != 0)
  {
    if (!lthread(t)) _kill(t->lt);
    if (!rthread(t)) _kill(t->rt);
    delete t;
  }
}


template <class SetType>
AVLSet<SetType>::AVLSet(AVLSet& b)
{
	cout << "AVLSet<SetType>::AVLSet(AVLSet& b)" << endl;
  if ((this->count = b.count) == 0)
  {
    root = 0;
  }
  else
  {
    _hold_nodes = new AVLSetNode<SetType>* [this->count];
    AVLSetNode<SetType>* t = b.leftmost();
    INT i = 0;
    while (t != 0)
    {
      _hold_nodes[i++] = new AVLSetNode<SetType>(t->item);
      t = b.succ(t);
    }
    root = _treeify(this->count);
  }
}


template <class SetType>
INT AVLSet<SetType>::operator == (AVLSet& y)
{
  if (this->count != y.count)
    return 0;
  else
  {
    AVLSetNode<SetType>* t = leftmost();
    AVLSetNode<SetType>* u = y.leftmost();
    for (;;)
    {
      if (t == 0)
        return 1;
      else if (!(SetTypeEQ(t->item, u->item)))
        return 0;
      else
      {
        t = succ(t);
        u = y.succ(u);
      }
    }
  }
}

template <class SetType>
INT AVLSet<SetType>::operator <= (AVLSet& y)
{
  if (this->count > y.count)
    return 0;
  else
  {
    AVLSetNode<SetType>* t = leftmost();
    AVLSetNode<SetType>* u = y.leftmost();
    for (;;)
    {
      if (t == 0)
        return 1;
      else if (u == 0)
        return 0;
      INT cmp = SetTypeCMP(t->item, u->item);
      if (cmp == 0)
      {
        t = succ(t);
        u = y.succ(u);
      }
      else if (cmp < 0)
        return 0;
      else
        u = y.succ(u);
    }
  }
}

template <class SetType>
void AVLSet<SetType>::operator |=(AVLSet& y)
{
  AVLSetNode<SetType>* t = leftmost();
  AVLSetNode<SetType>* u = y.leftmost();
  INT rsize = this->count + y.count;
  _hold_nodes = new AVLSetNode<SetType>* [rsize];
  INT k = 0;
  INT nzero=0;
  for (;;)
  {
    if (t == 0)
    {
      while (u != 0)
      {
//      	cout << "OOOOOOOOOOOOOOOOO:" << *u->item << endl;
        _hold_nodes[k++] = new AVLSetNode<SetType>(u->item, nzero);
//     	cout << "OOOOOOOOOOOOOOOOO:" << *_hold_nodes[k-1]->item << endl;
        u = y.succ(u);
      }
      break;
    }
    else if (u == 0)
    {
      while (t != 0)
      {
        _hold_nodes[k++] = t;
        t = succ(t);
      }
      break;
    }
    INT cmp = SetTypeCMP(t->item, u->item);
    if (cmp == 0)
    {
      _hold_nodes[k++] = t;
      t = succ(t);
      u = y.succ(u);
    }
    else if (cmp < 0)
    {
      _hold_nodes[k++] = t;
      t = succ(t);
    }
    else
    {
//      	cout << "UUUUUUUUUUUUUUUUU:" << *u->item << endl;
      _hold_nodes[k++] = new AVLSetNode<SetType>(u->item, nzero);
      u = y.succ(u);
    }
  }
  root = _treeify(k);
  this->count = k;
}

template <class SetType>
void AVLSet<SetType>::operator &= (AVLSet& y) 
{
  AVLSetNode<SetType>* t = leftmost();
  AVLSetNode<SetType>* u = y.leftmost();
//FD  INT rsize = (count < y.count)? count : y.count;
   INT rsize = y.count;
  _hold_nodes = new AVLSetNode<SetType>* [rsize];
  INT k = 0;
  for (;;)
  {
    if (t == 0)
      break;
    if (u == 0)
    {
      while (t != 0)
      {
        AVLSetNode<SetType>* tmp = succ(t);
        delete t;
        t = tmp;
      }
      break;
    }
    INT cmp = SetTypeCMP(t->item, u->item);
    if (cmp == 0)
    {
      _hold_nodes[k++] = t;
      t = succ(t);
      u = y.succ(u);
    }
    else if (cmp < 0)
    {
      AVLSetNode<SetType>* tmp = succ(t);
      delete t;
      t = tmp;
    }
    else
      u = y.succ(u);
  }
  root = _treeify(k);
  this->count = k;
}


template <class SetType>
void AVLSet<SetType>::operator -=(AVLSet& y)
{
  AVLSetNode<SetType>* t = leftmost();
  AVLSetNode<SetType>* u = y.leftmost();
  INT rsize = this->count;
  _hold_nodes = new AVLSetNode<SetType>* [rsize];
  INT k = 0;
  for (;;)
  {
    if (t == 0)
      break;
    else if (u == 0)
    {
      while (t != 0)
      {
        _hold_nodes[k++] = t;
        t = succ(t);
      }
      break;
    }
    INT cmp = SetTypeCMP(t->item, u->item);
    if (cmp == 0)
    {
      AVLSetNode<SetType>* tmp = succ(t);
      delete t;
      t = tmp;
      u = y.succ(u);
    }
    else if (cmp < 0)
    {
      _hold_nodes[k++] = t;
      t = succ(t);
    }
    else
      u = y.succ(u);
  }
  root = _treeify(k);
  this->count = k;
}

template <class SetType>
INT AVLSet<SetType>::owns(Pix i)
{
  if (i == 0) return 0;
  for (AVLSetNode<SetType>* t = leftmost(); t != 0; t = succ(t)) 
    if (Pix(t) == i) return 1;
  return 0;
}

template <class SetType>
INT AVLSet<SetType>::OK()
{
  INT v = 1;
  if (root == 0) 
    v = this->count == 0;
  else
  {
    INT n = 1;
    AVLSetNode<SetType>* trail = leftmost();
    AVLSetNode<SetType>* t = succ(trail);
    while (t != 0)
    {
      ++n;
      v &= SetTypeCMP(trail->item, t->item) < 0;
      trail = t;
      t = succ(t);
    }
    v &= n == this->count;
  }
  if (!v) this->error("invariant failure");
  return v;
}

/*
#include "MOType.h"

#include "Configuration.h"

#include "TupelStructureSel.h"


#include "InternalConfsSet.h"
#include "TupelStructureSet.h"
#include "extMOsSet.h"
#include "extEntry.h"
*/
/*
template <class SetType>
AVLSetNode<SetType>::AVLSetNode(SetType& h, INT dummy)
:lt(0), rt(0), stat(0)
{
	item = h->new_Set();
}
*/


#ifdef AVLSET_CLONE

template <class SetType>
AVLSetNode<SetType>::AVLSetNode(SetType& h, INT dummy)
//:item(h), lt(l), rt(r), stat(0)
:lt(0), rt(0), stat(0)
{
	item = (SetType) h->clone(0);
}
#else
template <class SetType>
AVLSetNode<SetType>::AVLSetNode(SetType& h, INT dummy)
:lt(0), rt(0), item(h), stat(0)
{
}
#endif
#undef AVLSET_CLONE

/*
template class AVLSet<Configuration<MOType> *>;

template class AVLSet<TupelStructureSel *>;

template class AVLSet<InternalConfsSet *>;
template class AVLSet<TupelStructureSet *>;
template class AVLSet<extMOsSet *>;
template class AVLSet<extEntry *>;

*/

//template class AVLSet<double>;

