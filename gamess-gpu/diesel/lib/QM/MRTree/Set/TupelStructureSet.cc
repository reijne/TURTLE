//***********************************************************************
//
//	Name:			TupelStructureSet.cc
//
//	Description:	
//
//	Author:			Michael Hanrath
//
//	Version:		0.1
//
//	Date:			30.05.1997
//
//
//
//
//***********************************************************************

#include "TupelStructureSet.h"
#include "extMOsSet.h"
#include "InternalConfsSet.h"
#include "NExternalsSet.h"

#include <stdio.h>


TupelStructureSet::TupelStructureSet(InternalConfsSet *parent) :
        Tree<InternalConfsSet>(parent),
	SetContainer<extMOsSet>(),
        TupelStructureBase()
{
}

TupelStructureSet::TupelStructureSet(istream &s) :
        Tree<InternalConfsSet>(NULL),
        SetContainer<extMOsSet>(s),
        TupelStructureBase(s)
{
	for ( ContainerIterator iter = first() ; !isLast(iter) ; next(iter) )
		operator [] (iter)->setParent(this);
}


TupelStructureSet::TupelStructureSet(IrRep irrep, 
		const Configuration<MOType> & conf,
		InternalConfsSet *parent, INT _ReferenceFlag) :
        Tree<InternalConfsSet>(parent),
	TupelStructureBase(
		irrep,
		conf,
		_ReferenceFlag)
{
}

TupelStructureSet::~TupelStructureSet()
{
	for ( ContainerIterator iter = first() ; !isLast(iter) ; next(iter) )
		if ( operator [] (iter) )
			delete operator [] (iter);
}


void	TupelStructureSet::add(const Configuration<MOType> &external)
{
	extMOsSet *t = new extMOsSet(
		external.getNumberOfOpenShells(), 
		external.getNumberOfClosedShells(), 
		this);
		
	if ( !contains(t) )
		SetContainerTree<InternalConfsSet, extMOsSet>::add(t);
	else
	{
		extMOsSet	*tt = t;
		t = operator() (seek(t));
		delete tt;
	}
	t->add(external);
}


void	TupelStructureSet::writeToStream(ostream & s) const
{
	Container<extMOsSet>::writeToStream(s);
	TupelStructureBase::writeToStream(s);
}



ostream& operator<<(ostream & s, const TupelStructureSet &base)
{
	s << ((TupelStructureBase &) base);
	return s;	
}


/*
void	TupelStructureSet::operator |= (TupelStructureSet& b)
{
	SetContainerTree<InternalConfsBase, extMOsSet>::operator |= (b);
	for ( ContainerIterator iter = b.first() ; !b.isLast(iter) ; b.next(iter) )
		if ( b[iter] )
			(*operator [] (seek(b[iter]))) |= (*b[iter]);
}
*/

#define AVLSET_CLONE
#include "../../../Container/AVLSet.cc"
template class AVLSet<TupelStructureSet *>;

#include "../../../Container/Set.cc"
template class Set<TupelStructureSet *>;
