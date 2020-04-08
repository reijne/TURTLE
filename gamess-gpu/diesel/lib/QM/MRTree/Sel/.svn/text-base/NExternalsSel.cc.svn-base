//***********************************************************************
//
//	Name:			NExternalsSel.cc
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

#include "NExternalsSel.h"
#include "InternalConfsSel.h"
#include "TupelStructureSel.h"
#include "extMOsSel.h"


#include "../../Configuration/Excitation.h"
#include "../../../Math/etc/BinomialCoefficient.h"

#include "../../IO/Fortran/Fort31File.h"

#include "../../MO/Iterators/MOList.h"
#include "../../MO/MRMOs.h"


#include <stdlib.h>
#include <fstream>
#include "../../../Container/String.h"


NExternalsSel::NExternalsSel(istream &s, const MOIrReps *moIrReps) :
        NExternalsBase<InternalConfsSel>(s, moIrReps),
        Tree<void>((void *) NULL),
        IndexedContainer<InternalConfsSel>(s)
{
	for ( ContainerIterator iter = first() ; !isLast(iter) ; next(iter) )
		operator [] (iter)->setParent(this);
}

NExternalsSel::NExternalsSel(
	const MRMOs *mrmos,
	INT NumberOfElectrons,
	INT Multiplicity,
	INT ExcitationLevel,
	INT NumberOfRoots,
	const ConfigurationSet &confSet) :
	NExternalsBase<InternalConfsSel>(
		mrmos, 
		NumberOfElectrons, 
		Multiplicity),
        Tree<void>((void *) NULL),
        IndexedContainer<InternalConfsSel>(ExcitationLevel+1)
{
	totalSymmetry = confSet(confSet.first()).calcIrRep(*mrmos);
	
InternalConfsSel	*intern = new InternalConfsSel (0, this);
	(*this)[0] = intern;
	for ( INT j=1 ; j<ExcitationLevel+1 ; j++ )
		(*this)[j] = NULL;
		
	RefSAFs = 0;
Pix	j = confSet.first();
	while ( j )
	{
	TupelStructureSel *tupel = new TupelStructureSel(
			0, totalSymmetry,
			confSet(j),
			intern,
			1);
		intern->add(tupel);
/*	 extMOsSel	*extMOs = new extMOsSel(1, 0, 0, tupel);
		(*tupel)[0] = extMOs;
		extMOs->setEnergy(RootEnergies(NumberOfRoots, 999));
*/
		(*tupel)[0]->setEnergy(RootEnergies(NumberOfRoots, 999));
		RefSAFs += (*eigFuncs)(confSet(j).getNumberOfOpenShells());
		confSet.next(j);
	}
	SAFs = init(0, mrmos->getMOSymmetryP());
}

NExternalsSel::~NExternalsSel()
{
	for ( ContainerIterator iter = first() ; !isLast(iter) ; next(iter) )
		if ( operator [] (iter) )
			delete operator [] (iter);
}

void	NExternalsSel::addExtEntry(
	const Configuration<MOType> & conf, const RootEnergies & r)
{
	(*this)[conf.getNumberOfExternals(getMRMOs())]->addExtEntry(conf, r);
}

const extEntry *	NExternalsSel::getExtEntry(
	Configuration<MOType> conf) const
{
//	cout << "HALLO: NExternalsSel::getExtEntry" << conf.getNumberOfExternals(getMRMOs()) << endl;
	return (*this)[conf.getNumberOfExternals(getMRMOs())]
		->getExtEntry(conf);
}

void	NExternalsSel::cutTreeLevel2()
{
	for ( ContainerIterator iter = first() ; !isLast(iter) ; next(iter) )
		if ( operator [] (iter) )
			operator [] (iter)->cutTreeLevel2();	
}

void	NExternalsSel::threshCut(EnergyType Threshold, INT divideByCSFs)
{
	for ( ContainerIterator iter = first() ; !isLast(iter) ; next(iter) )
		if ( operator [] (iter) )
			operator [] (iter)->threshCut(Threshold, divideByCSFs);	
}

void	NExternalsSel::writeToStream(ostream & s) const
{
	NExternalsBase<InternalConfsSel>::writeToStream(s);
	Container<InternalConfsSel>::writeToStream(s);
}

ostream& operator<<(ostream & s, const NExternalsSel &base)
{
	s << ((NExternalsBase<InternalConfsSel> &) base);
	return s;	
}

/*
istream& operator>>(istream & s, NExternalsSel &base)
{
	s >> ((NExternalsBase<InternalConfsSel> &) base);

INT	n;
	s >> n;
	for ( ContainerIterator iter = base.first() ; !base.isLast(iter) ;
		base.next(iter) )
	{
		base.operator[](iter) = new InternalConfsSel(&base);
		s >> *base.operator[](iter);
	}
	return s;	
}
*/
