//***********************************************************************
//
//	Name:			Container.h
//
//	Description:	abstract container
//
//	Author:			Michael Hanrath
//
//	Version:		0.0
//
//	Date:			14.03.1997
//
//
//
//
//
//***********************************************************************

#ifndef __Container_h
#define __Container_h

#include "../../../../config.h"

#include <iostream>
using std::ostream;
using std::istream;

union ContainerIterator;

template <class ContainedObjectType> class Container;
template <class ContainedObjectType> ostream& operator<< (ostream & s, const Container<ContainedObjectType> &);
template <class ContainedObjectType> istream& operator>> (istream & s, Container<ContainedObjectType> &);


template <class ContainedObjectType>
class Container {
public:
	Container() {}
//	Container(istream &s) {}
	
	
	virtual ContainerIterator	first() const = 0;
	virtual void	next(ContainerIterator &) const = 0;
	virtual INT	isLast(ContainerIterator) const = 0;
	
	virtual ContainedObjectType *& operator [] (ContainerIterator) = 0;
	virtual ContainedObjectType * const & operator [] (ContainerIterator) const = 0;

	virtual INT getNumberOfElements() const = 0;
	INT getNumberOfNotNILElements() const;

//--------------------------------------------------------------------------

	void	writeToStream(ostream & s) const;

	friend ostream& operator<< <ContainedObjectType> (ostream & s, const Container<ContainedObjectType> &);
	friend istream& operator>> <ContainedObjectType> (istream & s, Container<ContainedObjectType> &);

//--------------------------------------------------------------------------

private:
//	virtual Container<ContainedObjectType *> * new_Container();
};

template <class ContainedObjectType>	ostream& operator<<
	(ostream & s, const Container<ContainedObjectType> &);
template <class ContainedObjectType>	istream& operator>>
	(istream & s, Container<ContainedObjectType> &);


#endif
