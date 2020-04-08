//***********************************************************************
//
//	Name:			InternalConfsSet.cc
//
//	Description:	
//
//	Author:			Michael Hanrath
//
//	Version:		0.0
//
//	Date:			30.05.1997
//
//
//
//
//
//***********************************************************************

#include "InternalConfsSet.h"
#include "NExternalsSet.h"
#include "TupelStructureSet.h"
#include "extMOsSet.h"

#include <stdio.h>


/*
InternalConfsSet::InternalConfsSet(NExternalsSet *parent) :
	Tree<NExternalsSet>(parent)
{
	generated = 0;
}
*/

InternalConfsSet::InternalConfsSet(istream &s) :
        Tree<NExternalsSet>(NULL),
	SetContainer<TupelStructureSet>(s),
	InternalConfsBase(s)
{
	for ( ContainerIterator iter = first() ; !isLast(iter) ; next(iter) )
		operator [] (iter)->setParent(this);
	generated = 0;
}


InternalConfsSet::~InternalConfsSet()
{
	for ( ContainerIterator iter = first() ; !isLast(iter) ; next(iter) )
		delete operator [] (iter);
}


INT InternalConfsSet::checkSameRefs(const InternalConfsSet & a) const
{
	for ( ContainerIterator iter = first() ; !isLast(iter) ; next(iter) )
		if ( (*this)[iter]->isReference() )
		{
		TupelStructureSet *t = (TupelStructureSet *) (*this)[iter];
			if ( !a.seek(t) )
				return 0;
		}

	for ( ContainerIterator iter = a.first() ; !a.isLast(iter) ; a.next(iter) )
		if ( a[iter]->isReference() )
		{
		TupelStructureSet *t = (TupelStructureSet *) a[iter];
			if ( !seek(t) )
				return 0;
		}
		
	return 1;
}


void	InternalConfsSet::cutTreeLevel2()
{
	for ( ContainerIterator iter = first() ; !isLast(iter) ; next(iter) )
	{
		if ( operator [] (iter) )
		{
			if ( operator [] (iter)->getNumberOfLeaves()==0 )
			{
				delete operator [] (iter);
				operator [] (iter) = NULL;
			}
		}
	}
}


void	InternalConfsSet::add(const Configuration<MOType> &internal,
	const Configuration<MOType> &external)
{
	TupelStructureSet *t = new TupelStructureSet(
		internal.calcIrRep(*getParent()->getMRMOs()),
		internal, this);
		
	if ( !contains(t) )
		SetContainerTree<NExternalsSet, TupelStructureSet>::add(t);
	else
	{
		TupelStructureSet	*tt = t;
		t = operator() (seek(t));
		delete tt;
	}
	t->add(external);
}




/*
void	InternalConfsSet::setNumberOfConfigurations(INT n)
{
	number = n;
	mainConfiguration = new TupelStructureSet[number];
}
*/

void	InternalConfsSet::writeToStream(ostream & s) const
{
	Container<TupelStructureSet>::writeToStream(s);
	InternalConfsBase::writeToStream(s);
}

ostream& operator<<(ostream & s, const InternalConfsSet &base)
{
	s << ((InternalConfsBase &) base);
	s << ((Container<TupelStructureSet> &) base);
	return s;	
}

/*
istream& operator>>(istream & s, InternalConfsSet &base)
{
	s >> ((InternalConfsBase &) base);
	s >> ((Container<TupelStructureSet> &) base);
	return s;	
}
*/

/*
void	InternalConfsSet::operator |= (InternalConfsSet& b)
{
	SetContainerTree<NExternalsSet, TupelStructureSet>::operator |= (b);
	for ( ContainerIterator iter = b.first() ; !b.isLast(iter) ; b.next(iter) )
		if ( b[iter] )
			(*operator [] (seek(b[iter]))) |= (*b[iter]);
}
*/

#define AVLSET_CLONE
#include "../../../Container/AVLSet.cc"
template class AVLSet<InternalConfsSet *>;

#include "../../../Container/Set.cc"
template class Set<InternalConfsSet *>;
