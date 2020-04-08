//***********************************************************************
//
//	Name:			TupelStructureDiag.cc
//
//	Description:	
//
//	Author:			Michael Hanrath
//
//	Version:		0.1
//
//	Date:			21.06.1996
//					09.03.1997
//
//
//
//
//
//***********************************************************************

#include "TupelStructureDiag.h"
#include "extMOsDiag.h"
#include "InternalConfsDiag.h"
#include "NExternalsDiag.h"
#include "../Set/TupelStructureSet.h"
#include "../Set/extMOsSet.h"

#include <stdio.h>


TupelStructureDiag::TupelStructureDiag(const TupelStructureDiag &t) :
        Tree<InternalConfsDiag>(NULL),
	IndexedContainer<extMOsDiag>(t.getNumberOfElements()),
        TupelStructureBase(t)
{
	SAFStart = t.SAFStart;
	SAFInc = t.SAFInc;

	for ( ContainerIterator iter = t.first() , iter1 = first() ; 
		!t.isLast(iter) ; 
		t.next(iter), next(iter1) )
	{
		operator [] (iter1) = new extMOsDiag(*t[iter]);
		operator [] (iter1)->setParent(this);
	}
}


TupelStructureDiag::TupelStructureDiag(istream &s) :
        Tree<InternalConfsDiag>(NULL),
	IndexedContainer<extMOsDiag>(s),
        TupelStructureBase(s)
{
	for ( ContainerIterator iter = first() ; !isLast(iter) ; next(iter) )
		operator [] (iter)->setParent(this);
}

TupelStructureDiag::TupelStructureDiag(
	const TupelStructureSet &set, InternalConfsDiag *parent, INT referenceFlag) :
        Tree<InternalConfsDiag>(parent),
//	IndexedContainer<extMOsDiag>(set.length()),
	IndexedContainer<extMOsDiag>(parent->getNExt()/2+1),
        TupelStructureBase(parent->getParent()->getTotalSymmetry(),
                            (Configuration<MOType> &) set)
{
ContainerIterator j=set.first();
	for ( INT i=0 ; i<parent->getNExt()/2+1 ; i++ )
		p[i] = NULL;
	for ( INT i=0 ; 
		i<set.length() ; i++ , set.next(j) )
//		if ( parent->getNExt() )
		{
//			cout << "LLLLLLLLLL i=" << i << "-----" << *set[j] << endl;
//			cout << parent->getNExt() << " " << set[j]->getNumberOfOpenMOs()/2 << endl;
			p[set[j]->getNumberOfOpenMOs()/2] = new extMOsDiag(*set[j], this);
			ReferenceFlag = referenceFlag;
		}
}	

INT	TupelStructureDiag::init(INT n, const IrRep *irrep)
{
	SAFStart = n;
	SAFInc = getParent()->getParent()->getNumberOfSpinAdaptedFunctions(
			getNumberOfOpenShells());

	if ( getNumberOfElements() )
		return MRTreeBase<InternalConfsDiag, extMOsDiag>::
			init(n, irrep);

	return 0;
}

/*
void	TupelStructureDiag::setNumber(INT n)
{
	if ( number == n )
		return;
	delete moContainer;
	number = n;
}
*/


void	TupelStructureDiag::writeToStream(ostream & s) const
{
	Container<extMOsDiag>::writeToStream(s);
	TupelStructureBase::writeToStream(s);
}



ostream& operator<<(ostream & s, const TupelStructureDiag &base)
{
	s << ((TupelStructureBase &) base);
	s << ((Container<extMOsDiag> &) base);
	return s;	
}

/*
istream& operator>>(istream & s, TupelStructureDiag &base)
{
	s >> ((TupelStructureBase &) base);
	s >> ((Container<extMOsDiag> &) base);
	return s;	
}
*/
