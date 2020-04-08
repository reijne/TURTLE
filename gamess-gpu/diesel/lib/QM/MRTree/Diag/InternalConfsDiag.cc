//***********************************************************************
//
//	Name:			InternalConfsDiag.cc
//
//	Description:	
//
//	Author:			Michael Hanrath
//
//	Version:		0.0
//
//	Date:			21.06.1996
//
//
//
//
//
//***********************************************************************

#include "InternalConfsDiag.h"
#include "TupelStructureDiag.h"
#include "../Set/InternalConfsSet.h"

#include <stdio.h>


InternalConfsDiag::InternalConfsDiag(const InternalConfsDiag &t) :
        Tree<NExternalsDiag>(NULL),
	IndexedContainer<TupelStructureDiag>(t.getNumberOfElements()),
        InternalConfsBase(t)
{
	for ( ContainerIterator iter = t.first() , iter1 = first() ; 
		!t.isLast(iter) ; 
		t.next(iter), next(iter1) )
	{
		operator [] (iter1) = new TupelStructureDiag(*t[iter]);
		operator [] (iter1)->setParent(this);
	}
}


InternalConfsDiag::InternalConfsDiag(istream &s) :
        Tree<NExternalsDiag>(NULL),
	IndexedContainer<TupelStructureDiag>(s),
        InternalConfsBase(s)
{
	for ( ContainerIterator iter = first() ; !isLast(iter) ; next(iter) )
		operator [] (iter)->setParent(this);
}

/*
void	InternalConfsDiag::setNumberOfConfigurations(INT n)
{
	number = n;
	mainConfiguration = new TupelStructureDiag[number];
}
*/

InternalConfsDiag::InternalConfsDiag(
	const InternalConfsSet &set, NExternalsDiag *parent, INT referenceFlag) :
        Tree<NExternalsDiag>(parent),
	IndexedContainer<TupelStructureDiag>(set.length()),
        InternalConfsBase(set.getNExt())
{
ContainerIterator j=set.first();
	for ( INT i=0 ; 
		i<getNumberOfElements() ; i++ , set.next(j) )
	{
//		cout << "i= " << i << endl;
		p[i] = new TupelStructureDiag(*set[j], this, referenceFlag);
	}
}
	


void	InternalConfsDiag::writeToStream(ostream & s) const
{
	Container<TupelStructureDiag>::writeToStream(s);
	InternalConfsBase::writeToStream(s);
}

ostream& operator<<(ostream & s, const InternalConfsDiag &base)
{
	s << ((InternalConfsBase &) base);
	return s;	
}

/*
istream& operator>>(istream & s, InternalConfsDiag &base)
{
	s >> ((InternalConfsBase &) base);
	s >> ((Container<TupelStructureDiag> &) base);
	return s;	
}
*/
