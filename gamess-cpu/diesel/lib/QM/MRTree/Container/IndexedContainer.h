//***********************************************************************
//
//	Name:			IndexedContainer.h
//
//	Description:	
//
//	Author:			Michael Hanrath
//
//	Version:		0.0
//
//	Date:			30. Aug 1998
//
//***********************************************************************

#ifndef __INDEXEDCONTAINER_H
#define __INDEXEDCONTAINER_H

#include "../../../../config.h"

#include <iostream>
using std::cout;

#include "Container.h"
#include "ContainerIterator.h"

template <class ContainedObjectType>
class IndexedContainer;

template <class ContainedObjectType>
ostream& operator<< (ostream & s, const IndexedContainer<ContainedObjectType> &);

template <class ContainedObjectType>
istream& operator>> (istream & s, IndexedContainer<ContainedObjectType> &);

template <class ContainedObjectType>
class IndexedContainer :
	virtual public Container<ContainedObjectType> {
public:
	IndexedContainer() { cout << "Error: default constructor IndexedContainer::IndexedContainer() called"; exit(1); }
	IndexedContainer(INT n);
	IndexedContainer(istream &s);
	virtual ~IndexedContainer();
	
	ContainerIterator	first() const;
	void	next(ContainerIterator &) const;
	INT	isLast(ContainerIterator) const;
	
	ContainedObjectType *& operator [] (ContainerIterator);
	ContainedObjectType *& operator [] (INT i);

	ContainedObjectType * const & operator [] (ContainerIterator) const;
	ContainedObjectType * const & operator [] (INT i) const;

	INT getNumberOfElements() const;


//--------------------------------------------------------------------------

	friend ostream& operator<< <> (ostream & s, const IndexedContainer<ContainedObjectType> &);
	friend istream& operator>> <> (istream & s, IndexedContainer<ContainedObjectType> &);

//--------------------------------------------------------------------------


protected:
ContainedObjectType	**p;
INT	n;
};

/*
template <class ContainedObjectType>	ostream& operator<<
	(ostream & s, const IndexedContainer<ContainedObjectType> &);
template <class ContainedObjectType>	istream& operator>>
	(istream & s, IndexedContainer<ContainedObjectType> &);
*/


template <class ContainedObjectType>
inline 
ContainerIterator	IndexedContainer<ContainedObjectType>::first() const
{INT nzero=0;	return	nzero;	}

template <class ContainedObjectType>
inline 
void	IndexedContainer<ContainedObjectType>::next(ContainerIterator &iterator) const
{	iterator.i++;	}

template <class ContainedObjectType>
inline 
INT	IndexedContainer<ContainedObjectType>::isLast(ContainerIterator iterator) const
{	return (iterator.i>=n);	}

template <class ContainedObjectType>
inline 
ContainedObjectType *& IndexedContainer<ContainedObjectType>
	::operator [] (ContainerIterator iterator)
{	return p[iterator.i];	}


template <class ContainedObjectType>
inline ContainedObjectType *& IndexedContainer<ContainedObjectType>
	::operator [] (INT i)
{	return p[i];	}

template <class ContainedObjectType>
inline 
ContainedObjectType * const & IndexedContainer<ContainedObjectType>
	::operator [] (ContainerIterator iterator) const
{	return p[iterator.i];	}


template <class ContainedObjectType>
inline 
ContainedObjectType * const & IndexedContainer<ContainedObjectType>
	::operator [] (INT i) const
{	return p[i];	}

template <class ContainedObjectType>
inline INT IndexedContainer<ContainedObjectType>::getNumberOfElements() const
{	return n;	}


#endif
