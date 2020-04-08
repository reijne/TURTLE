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


#ifndef _AVLMap_h
#define _AVLMap_h 1

#include "../../config.h"

#include "Map.h"

template <class KeyType, class ObjectType>
struct AVLNode
{
  AVLNode<KeyType, ObjectType>*      lt;
  AVLNode<KeyType, ObjectType>*      rt;
  KeyType                 item;
  ObjectType                 cont;
  char                stat;
                      AVLNode(KeyType  h, ObjectType  c, 
                                    AVLNode<KeyType, ObjectType>* l=0, AVLNode<KeyType, ObjectType>* r=0);
                      ~AVLNode();
};

template <class KeyType, class ObjectType>
inline AVLNode<KeyType, ObjectType>::AVLNode(KeyType  h, ObjectType  c, 
                                    AVLNode<KeyType, ObjectType>* l, AVLNode<KeyType, ObjectType>* r)
     :lt(l), rt(r), item(h), cont(c), stat(0) {}

template <class KeyType, class ObjectType>
inline AVLNode<KeyType, ObjectType>::~AVLNode() {}

//typedef AVLNode<KeyType, ObjectType>* AVLNodePtr;


template <class KeyType, class ObjectType>
class AVLMap : public Map<KeyType, ObjectType>
{
private:
static KeyType*   _target_item;     // add/del_item target
static AVLNode<KeyType, ObjectType>* _found_node; // returned added/deleted node
protected:
  AVLNode<KeyType, ObjectType>*   root;

  AVLNode<KeyType, ObjectType>*   leftmost();
  AVLNode<KeyType, ObjectType>*   rightmost();
  AVLNode<KeyType, ObjectType>*   pred(AVLNode<KeyType, ObjectType>* t);
  AVLNode<KeyType, ObjectType>*   succ(AVLNode<KeyType, ObjectType>* t);
  void            _kill(AVLNode<KeyType, ObjectType>* t);
  void            _add(AVLNode<KeyType, ObjectType>*& t);
  void            _del(AVLNode<KeyType, ObjectType>* p, AVLNode<KeyType, ObjectType>*& t);

public:
                AVLMap(ObjectType  dflt);
                AVLMap(AVLMap<KeyType, ObjectType>& a);
                ~AVLMap();

  ObjectType&          operator [] (KeyType  key);

  void          del(KeyType  key);

  Pix           first();
  void          next(Pix& i);
  KeyType&          key(Pix i);
  ObjectType&          contents(Pix i);

  Pix           seek(KeyType  key);
  INT           contains(KeyType  key);

  void          clear(); 

  Pix           last();
  void          prev(Pix& i);

  INT           OK();
};

template <class KeyType, class ObjectType>
inline AVLMap<KeyType, ObjectType>::~AVLMap()
{
  _kill(root);
}

template <class KeyType, class ObjectType>
inline AVLMap<KeyType, ObjectType>::AVLMap(ObjectType  dflt) :Map<KeyType, ObjectType>(dflt)
{
  root = 0;
}

template <class KeyType, class ObjectType>
inline Pix AVLMap<KeyType, ObjectType>::first()
{
  return Pix(leftmost());
}

template <class KeyType, class ObjectType>
inline Pix AVLMap<KeyType, ObjectType>::last()
{
  return Pix(rightmost());
}

template <class KeyType, class ObjectType>
inline void AVLMap<KeyType, ObjectType>::next(Pix& i)
{
  if (i != 0) i = Pix(succ((AVLNode<KeyType, ObjectType>*)i));
}

template <class KeyType, class ObjectType>
inline void AVLMap<KeyType, ObjectType>::prev(Pix& i)
{
  if (i != 0) i = Pix(pred((AVLNode<KeyType, ObjectType>*)i));
}

template <class KeyType, class ObjectType>
inline KeyType& AVLMap<KeyType, ObjectType>::key(Pix i)
{
  if (i == 0) this->error("null Pix");
  return ((AVLNode<KeyType, ObjectType>*)i)->item;
}

template <class KeyType, class ObjectType>
inline ObjectType& AVLMap<KeyType, ObjectType>::contents(Pix i)
{
  if (i == 0) this->error("null Pix");
  return ((AVLNode<KeyType, ObjectType>*)i)->cont;
}

template <class KeyType, class ObjectType>
inline void AVLMap<KeyType, ObjectType>::clear()
{
  _kill(root);
  this->count = 0;
  root = 0;
}

template <class KeyType, class ObjectType>
inline INT AVLMap<KeyType, ObjectType>::contains(KeyType  key)
{
  return seek(key) != 0;
}



#endif
