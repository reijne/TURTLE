//***********************************************************************
//
//	Name:			NExternalsSet.cc
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

#include "NExternalsSet.h"
#include "InternalConfsSet.h"
#include "TupelStructureSet.h"
#include "extMOsSet.h"


#include "../../Configuration/Excitation.h"
#include "../../../Math/etc/BinomialCoefficient.h"

#include "../../IO/Fortran/Fort31File.h"

#include "../../MO/Iterators/MOList.h"


#include <stdlib.h>
#include <fstream>
#include "../../../Container/String.h"


NExternalsSet::NExternalsSet() :
        NExternalsBase<InternalConfsSet>(),
        Tree<void>((void *) NULL),
        SetContainer<InternalConfsSet>()
{
}

NExternalsSet::NExternalsSet(istream &s, const MRMOs *mrmos) :
        NExternalsBase<InternalConfsSet>(s, mrmos),
        Tree<void>((void *) NULL),
        SetContainer<InternalConfsSet>(s)
{
	for ( ContainerIterator iter = first() ; !isLast(iter) ; next(iter) )
		operator [] (iter)->setParent(this);
}

NExternalsSet::NExternalsSet(
	const MRMOs *mrmos,
	INT NumberOfElectrons,
	INT Multiplicity,
	INT NumberOfRoots,
	const ConfigurationSet &refConfs) :
	NExternalsBase<InternalConfsSet>(
		mrmos, 
		NumberOfElectrons, 
		Multiplicity), 
	Tree<void>((void *) NULL)
{
	totalSymmetry = refConfs(refConfs.first()).calcIrRep(*mrmos);
	
	setNumberOfRoots(NumberOfRoots);
	for ( INT i=0 ; i<NumberOfRoots ; i++ )
		setRootNumber(i, i+1);
	RefSAFs = 0;

InternalConfsSet	*intern = new InternalConfsSet (0, this);
	add(intern);

	Pix	i = refConfs.first();
	while ( i )
	{
	TupelStructureSet *tupel = new TupelStructureSet(
			totalSymmetry,
			refConfs(i),
			intern,
			1);
		intern->add(tupel);
		tupel->add(Configuration<MOType>());


		RefSAFs += (*eigFuncs)(refConfs(i).getNumberOfOpenShells());
		refConfs.next(i);
	}
	SAFs = init(0, mrmos->getMOSymmetryP());
	References = refConfs.length();
/*	
InternalConfsSet	*intern = new InternalConfsSet (0, this);
	(*this)[0] = intern;
	for ( INT i=1 ; i<ExcitationLevel+1 ; i++ )
		(*this)[i] = NULL;
		
	RefSAFs = 0;
Pix	i = refConfs.first();
	while ( i )
	{
	TupelStructureSel *tupel = new TupelStructureSel(
			0, totalSymmetry,
			*refConfs(i),
			intern,
			1);
		intern->add(tupel);

		(*tupel)[0]->setEnergy(RootEnergies(NumberOfRoots, 999));
		RefSAFs += (*eigFuncs)(refConfs(i)->getNumberOfOpenShells());
		refConfs.next(i);
	}
	SAFs = init(0, mrmos->getMOSymmetryP());
*/
//	initAll();
}




NExternalsSet::NExternalsSet(
	const MRMOs *mrmos,
	INT NumberOfElectrons,
	INT Multiplicity,
	INT NumberOfRoots,
	IrRep _totalSymmetry) :
	NExternalsBase<InternalConfsSet>(
		mrmos, 
		NumberOfElectrons, 
		Multiplicity), 
	Tree<void>((void *) NULL)
{
	totalSymmetry = _totalSymmetry;
	
	setNumberOfRoots(NumberOfRoots);
	setRootNumber(0, 1);
	RefSAFs = 0;

InternalConfsSet	*intern = new InternalConfsSet (0, this);
	add(intern);

	SAFs = init(0, mrmos->getMOSymmetryP());
}

NExternalsSet::~NExternalsSet()
{
	for ( ContainerIterator iter = first() ; !isLast(iter) ; next(iter) )
		if ( operator [] (iter) )
			delete operator [] (iter);
}


void	NExternalsSet::add(const Configuration<MOType> &conf)
{
	// check symmetry
	if ( length() )
	{
		if ( totalSymmetry != conf.calcIrRep(*mrmos) )
		{
                  std::cerr << "NExternalsSet::add: wrong symmetry." << std::endl;
			exit(1);
		}
	}
	else
		totalSymmetry = conf.calcIrRep(*mrmos);


	// split in internal/external part
Configuration<MOType>	internal;
Configuration<MOType>	external;

	conf.split(mrmos, internal, external);


	InternalConfsSet	*t = new InternalConfsSet(
		external.getNumberOfElectrons(), this);
		
	if ( !contains(t) )
		SetContainerTree<void, InternalConfsSet>::add(t);
	else
	{
		InternalConfsSet	*tt = t;
		t = operator() (seek(t));
		delete tt;
	}
	t->add(internal,external);

}


void	NExternalsSet::cutTreeLevel2()
{
	for ( ContainerIterator iter = first() ; !isLast(iter) ; next(iter) )
		if ( operator [] (iter) )
			operator [] (iter)->cutTreeLevel2();	
}

INT NExternalsSet::checkSameRefs(const NExternalsSet & ext) const
{
	if ( operator [] (first()) && ext[ext.first()] )
		return operator [] (first())->checkSameRefs(*ext[ext.first()]);
	return 0;
}


void	NExternalsSet::writeToStream(ostream & s) const
{
	NExternalsBase<InternalConfsSet>::writeToStream(s);
	Container<InternalConfsSet>::writeToStream(s);
}

ostream& operator<<(ostream & s, const NExternalsSet &base)
{
	s << ((NExternalsBase<InternalConfsSet> &) base);
	return s;	
}

/*
void	NExternalsSet::operator |= (NExternalsSet& b)
{
	SetContainerTree<void, InternalConfsSet>::operator |= (b);
	for ( ContainerIterator iter = b.first() ; !b.isLast(iter) ; b.next(iter) )
		if ( b[iter] )
			(*operator [] (seek(b[iter]))) |= (*b[iter]);
}
*/
