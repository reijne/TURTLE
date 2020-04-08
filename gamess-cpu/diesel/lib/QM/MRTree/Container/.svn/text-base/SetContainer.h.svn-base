//***********************************************************************
//
//	Name:			SetContainer.h
//
//	Description:	implements the functionality of a set based
//					on the GNU-C++-Library
//
//	Author:			Michael Hanrath
//
//	Version:		0.0
//
//	Date:			08.03.1997
//
//
//
//
//
//***********************************************************************


#ifndef __SETCONTAINER_H
#define __SETCONTAINER_H

#include "../../../../config.h"


#include "Container.h"
#include "ContainerIterator.h"
#include "../../../Container/AVLSet.h"


template <class SetType>
class _AVLSet {
public:

  INT                   length() const {return avl.length();}               // current number of items
  	
  INT                   empty() const{return avl.empty();}

  Pix           add(SetType& item) {return avl.add(item);}     // add item; return Pix
  void          del(SetType& item) {avl.del(item);}     // delete item
  INT           contains(SetType& item) {return avl.contains(item);}    // is item in set?

  void          clear() {avl.clear();}                // delete all items

  Pix           first() const{return avl.first();}              // Pix of first item or 0
  void          next(Pix& i) const {avl.next(i);}        // advance to next or 0
  SetType&          operator () (Pix i) const {return avl.operator() (i);} // access item at i

  INT           owns(Pix i) {return avl.owns(i);}             // is i a valid Pix  ?
  Pix           seek(SetType& item) const {return avl.seek(item);}         // Pix of item

  void                  operator |= (_AVLSet<SetType>& b) {avl |= b.avl;} // add all items in b
  void                  operator -= (_AVLSet<SetType>& b) {avl -= b.avl;} // delete items also in b
  void                  operator &= (_AVLSet<SetType>& b) {avl &= b.avl;} // delete items not in b

  INT                   operator == (_AVLSet<SetType>& b) {return avl==b.avl;};
  INT                   operator != (_AVLSet<SetType>& b) {return avl!=b.avl;}
  INT                   operator <= (_AVLSet<SetType>& b) {return avl<=b.avl;} 

  void                  error(const char* msg) const {avl.error(msg);}
  INT           OK() {return avl.OK();}                // rep invariant

private:
AVLSet<SetType>	avl;
};


template <class ContainedObjectType>
class SetContainer : 
	public virtual Container<ContainedObjectType>,
	public virtual _AVLSet<ContainedObjectType *> {
public:
	SetContainer() {}
//	SetContainer(INT n) {}
	SetContainer(istream &s);
	virtual ~SetContainer() {}

	// "virtual" constructor:
	virtual SetContainer<ContainedObjectType> * new_Set();
	
	ContainerIterator	first() const;
	void	next(ContainerIterator &) const;
	INT	isLast(ContainerIterator) const;
	
	ContainedObjectType *& operator [] (ContainerIterator);
	ContainedObjectType * const & operator [] (ContainerIterator) const;

	INT getNumberOfElements() const;

};


template <class ContainedObjectType>
inline 
SetContainer<ContainedObjectType> *	SetContainer<ContainedObjectType>::
	new_Set()
{	return new SetContainer<ContainedObjectType>();	}

template <class ContainedObjectType>
inline 
ContainerIterator	SetContainer<ContainedObjectType>::first() const
{	return	_AVLSet<ContainedObjectType *>::first();	}

template <class ContainedObjectType>
inline 
void	SetContainer<ContainedObjectType>::next(ContainerIterator &iterator) const
{	_AVLSet<ContainedObjectType *>::next(iterator.pix);	}

template <class ContainedObjectType>
inline 
INT	SetContainer<ContainedObjectType>::isLast(ContainerIterator iterator) const
{	return (iterator.pix == 0);	}

template <class ContainedObjectType>
inline 
ContainedObjectType *& SetContainer<ContainedObjectType>
	::operator [] (ContainerIterator iterator)
{	return _AVLSet<ContainedObjectType*>::operator () (iterator.pix);	}

template <class ContainedObjectType>
inline 
ContainedObjectType * const & SetContainer<ContainedObjectType>
	::operator [] (ContainerIterator iterator) const
{	return _AVLSet<ContainedObjectType *>::operator () (iterator.pix);	}

template <class ContainedObjectType>
inline 
INT SetContainer<ContainedObjectType>::getNumberOfElements() const
{	return this->length();	}




#endif
