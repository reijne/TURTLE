//***********************************************************************
//
//	Name:			TupelStructureSel.cc
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

#include "TupelStructureSel.h"
#include "extMOsSel.h"
#include "InternalConfsSel.h"
#include "NExternalsSel.h"

#include <stdio.h>

using namespace std;

TupelStructureSel::TupelStructureSel(InternalConfsSel *parent) :
        Tree<InternalConfsSel>(parent),
	IndexedContainer<extMOsSel>(),
        TupelStructureBase()
{
}

TupelStructureSel::TupelStructureSel(istream &s) :
        Tree<InternalConfsSel>(NULL),
	IndexedContainer<extMOsSel>(s),
        TupelStructureBase(s)
{
	for ( ContainerIterator iter = first() ; !isLast(iter) ; next(iter) )
		operator [] (iter)->setParent(this);
}

TupelStructureSel::TupelStructureSel(INT nElectrons, IrRep irrep, 
		Configuration<MOType> conf,
		InternalConfsSel *parent, INT _ReferenceFlag) :
        Tree<InternalConfsSel>(parent),
        IndexedContainer<extMOsSel>(nElectrons/2+1),
	TupelStructureBase(
		irrep,
		conf,
		_ReferenceFlag)
{
	for ( INT i=0 ; i<getNumberOfElements() ; i++ )
		(*this)[i] = 		
			new extMOsSel(1,
			2*i + (nElectrons & 1),
			nElectrons/2 - i,
			this);

}

TupelStructureSel::~TupelStructureSel()
{
	for ( ContainerIterator iter = first() ; !isLast(iter) ; next(iter) )
		if ( operator [] (iter) )
			delete operator [] (iter);
}


INT	TupelStructureSel::init(INT n, const IrRep *irrep)
{
	SAFStart = n;
	SAFInc = getParent()->getParent()->getNumberOfSpinAdaptedFunctions(
			getNumberOfOpenShells());

	if ( getNumberOfElements() )
		return MRTreeBase<InternalConfsSel, extMOsSel>::
			init(n, irrep);
	return 0;
}


/*
void	TupelStructureSel::add(
	const Configuration<MOType> & conf,
	INT	roots,
	const EnergyType *energy
	)
{
	if ( !this->operator [] (conf.getNumberOfOpenShells()/2) )
		this->operator [] (conf.getNumberOfOpenShells()/2) =
			new extMOsSel(1,
			conf.getNumberOfOpenShells(),
			conf.getNumberOfClosedShells(),
			this);
	this->operator [] (conf.getNumberOfOpenShells()/2)->add(
		conf, roots, energy);
}
*/

/*
void	TupelStructureSel::setNumber(INT n)
{
	if ( number == n )
		return;
	delete moContainer;
	number = n;
}
*/

void	TupelStructureSel::addExtEntry(
	const Configuration<MOType> & external, const RootEnergies & r)
{
//	cout << "HALLO: TupelStructureSel::addExtEntry" << external.getNumberOfOpenShells()/2 << endl;
	if ( external.getNumberOfElectrons()>n )
		return;
	(*this)[external.getNumberOfOpenShells()/2]->add(external, r);
	return;
}


const extEntry *	TupelStructureSel::getExtEntry(
	Configuration<MOType> external) const
{
//	cout << "HALLO: TupelStructureSel::getExtEntry" << external.getNumberOfOpenShells()/2 << endl;
	if ( external.getNumberOfElectrons()>n )
		return NULL;
	return (*this)[external.getNumberOfOpenShells()/2]->getExtEntry(external);
}



void	TupelStructureSel::threshCut(EnergyType Threshold, INT divideByCSFs)
{
	for ( ContainerIterator iter = first() ; !isLast(iter) ; next(iter) )
		if ( operator [] (iter) )
			operator [] (iter)->threshCut(Threshold, divideByCSFs);	
}


void	TupelStructureSel::writeToStream(ostream & s) const
{
	Container<extMOsSel>::writeToStream(s);
	TupelStructureBase::writeToStream(s);
}



ostream& operator<<(ostream & s, const TupelStructureSel &base)
{
	s << ((TupelStructureBase &) base);
	return s;	
}

/*
istream& operator>>(istream & s, TupelStructureSel &base)
{
	s >> ((TupelStructureBase &) base);
	s >> ((Container<extMOsSel> &) base);
	return s;	
}
*/

#define AVLSET_CLONE
#include "../../../Container/AVLSet.cc"
template class AVLSet<TupelStructureSel *>;

#include "../../../Container/Set.cc"
template class Set<TupelStructureSel *>;
