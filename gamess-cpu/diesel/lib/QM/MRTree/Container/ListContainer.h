//***********************************************************************
//
//	Name:			ListContainer.h
//
//	Description:	implements the functionality of a singly linked list
//					based on the GNU-C++-Library
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


#ifndef __ListContainer_h
#define __ListContainer_h

#include "../../../../config.h"


#include "Container.h"
#include "ContainerIterator.h"
#include "../../../Container/SLList.h"

template <class ContainedObjectType>
class ListContainer : 
	public virtual Container<ContainedObjectType>,
	public virtual SLList<ContainedObjectType *> {
public:
	ListContainer() {}
//	ListContainer(INT n) {}
	ListContainer(istream &s) :
		SLList<ContainedObjectType *> ()
{
//	cout << "ListContainer(istream &s)" << endl;
INT	n;
	s >> n;
//	cout << "n=" << n << endl;
	
	for ( INT i=0 ; i<n ; i++ )
	{
	ContainedObjectType	*p = new ContainedObjectType(s);
		append(p);
	}
//	cout << "OK1" << endl;
}


	ContainerIterator	first() const;
	void	next(ContainerIterator &) const;
	INT	isLast(ContainerIterator) const;
	
	ContainedObjectType *& operator [] (ContainerIterator);
	ContainedObjectType * const & operator [] (ContainerIterator) const;

	INT getNumberOfElements() const;
};

template <class ContainedObjectType>
inline 
ContainerIterator	ListContainer<ContainedObjectType>::first() const
{	return	SLList<ContainedObjectType *>::first();	}

template <class ContainedObjectType>
inline 
void	ListContainer<ContainedObjectType>::next(ContainerIterator &iterator) const
{	SLList<ContainedObjectType *>::next(iterator.pix);	}

template <class ContainedObjectType>
inline 
INT	ListContainer<ContainedObjectType>::isLast(ContainerIterator iterator) const
{	return (iterator.pix == 0);	}

template <class ContainedObjectType>
inline 
ContainedObjectType *& ListContainer<ContainedObjectType>
	::operator [] (ContainerIterator iterator)
{	return SLList<ContainedObjectType *>::operator () (iterator.pix);	}

template <class ContainedObjectType>
inline 
ContainedObjectType * const & ListContainer<ContainedObjectType>
	::operator [] (ContainerIterator iterator) const
{	return SLList<ContainedObjectType *>::operator () (iterator.pix);	}

template <class ContainedObjectType>
inline 
INT ListContainer<ContainedObjectType>::getNumberOfElements() const
{	return this->length();	}




#endif
