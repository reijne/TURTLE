//***********************************************************************
//
//	Name:			extMOsSet.cc
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

#include <stdlib.h>
#include <string>






extMOsSet::~extMOsSet()
{
}

void	extMOsSet::add(const Configuration<MOType> &external)
{
Configuration<MOType>	conf;
	conf = external;
	extEntry *t = new extEntry(this, 0, NULL, 0, NULL, conf, 0);
		
	if ( !contains(t) )
		SetContainerTree<TupelStructureBase, extEntry>::add(t);
	else
	{
		extEntry	*tt = t;
		t = operator() (seek(t));
		delete tt;
	}
}



INT	extMOsSet::init(INT n, const IrRep *irrep)
{
	// initialize SAFStart
	SAFStart = n;
	SAFInc = getParent()->getParent()->getParent()->
			getNumberOfSpinAdaptedFunctions(
			open + getParent()->getNumberOfOpenShells());
			
	return n + length()*SAFInc;
}


void	extMOsSet::writeToStream(ostream & s) const
{
//	s << "Container<extEntry>::writeToStream(s)" << endl;
	Container<extEntry>::writeToStream(s);
//	s << "extMOsBase::writeToStream(s)" << endl;
	extMOsBase::writeToStream(s);
}

ostream& operator<<(ostream & s, const extMOsSet & base)
{
	s << ((extMOsBase &) base);
	return s;	


}


ConfigurationSAFNr<MOType>	extMOsSet::
	getConfigurationSAFNr(INT ind) const
{
ContainerIterator	iter = first();

	for ( INT i=1 ; i<ind && !isLast(iter) ; i++ )
		next(iter);

	return getConfigurationSAFNr(iter.pix);
}


ConfigurationSAFNr<MOType>	extMOsSet::getConfigurationSAFNr(Pix pix) const
{
Configuration<MOType>	conf = *getParent();


//	cout << "n= " << n << ", this = " << this << endl;
	for ( INT i=0 ; i<getNumberOfOpenMOs() ; i++ )
		conf.create(getOpenMO(pix, i));

	for ( INT i=0 ; i<getNumberOfClosedMOs() ; i++ )
	{	conf.create(getClosedMO(pix, i));
		conf.create(getClosedMO(pix, i));
	}


INT	j = 0;

Pix	ii = _AVLSet<extEntry *>::first();
	while ( ii )
	{
		if ( ii==pix )
			break;
		_AVLSet<extEntry *>::next(ii);
		j++;
	}


	return ConfigurationSAFNr<MOType>(
		conf, *getParent(),
		getNumberOfOpenMOs(), getNumberOfClosedMOs(),
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		SAFStart + j*SAFInc, SAFInc,
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		getEnergy(pix),
		getParent()->isReference(),
		getSelectionFlags(pix)
		);

}


/*
void	extMOsSet::operator |= (extMOsSet& b)
{
	SetContainerTree<TupelStructureBase, extEntry>::operator |= (b);
}
*/

#define AVLSET_CLONE
#include "../../../Container/AVLSet.cc"
template class AVLSet<extMOsSet *>;

#include "../../../Container/Set.cc"
template class Set<extMOsSet *>;
