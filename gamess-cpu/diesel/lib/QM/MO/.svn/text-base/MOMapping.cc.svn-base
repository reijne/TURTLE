//***********************************************************************
//
//	Name:			MOMapping.cc
//
//	Description:	
//
//	Author:			Michael Hanrath
//
//	Version:		0.0
//
//	Date:			12.05.1997
//
//
//
//
//
//***********************************************************************

#include "MOMapping.h"

#include <iomanip>


#include "../MRTree/Diag/NExternalsDiag.h"
#include "../MRTree/Set/NExternalsSet.h"
#include "../RepresentationMatrices/CIVectors.h"

#include "../MRTree/MRTreeIterator.h"
#include "../MRTree/Set/extMOsSet.h"

#include "../../Container/AVLMap.h"

using namespace std;

MOMapping	moMapping;

MOMapping::MOMapping()
{
	maxMO = 0;
	mapRealToContinuous = NULL;
	mapContinuousToReal = NULL;
}


MOMapping::MOMapping(istream &s)
{
	s >> maxMO;
	mapRealToContinuous = new MOType[maxMO];
	mapContinuousToReal = new MOType[maxMO];
	for ( INT i=0 ; i<maxMO ; i++ )
		s >> mapRealToContinuous[i];
	for ( INT i=0 ; i<maxMO ; i++ )
		s >> mapContinuousToReal[i];
}


MOMapping::MOMapping(INT _maxMO, INT *internal, const MOIrReps &moirreps)
{
	maxMO = _maxMO;
	mapRealToContinuous = new MOType[maxMO];
	mapContinuousToReal = new MOType[maxMO];
INT	n = 0;
	for ( IrRep irrep=0 ; irrep<moirreps.getNumberOfIrReps() ; irrep++ )
	{
	MOType mo = moirreps.getStartMO(irrep);
	MOType moi = mo;
		for ( INT i=0 ; i<moirreps.getInIrRep(irrep) ; i++ , moi++ )
			if ( internal[n+i] )
			{
				mapRealToContinuous[moi-1] = mo;
				mapContinuousToReal[mo-1] = moi;
				mo++;
			}
		moi = moirreps.getStartMO(irrep);
		for ( INT i=0 ; i<moirreps.getInIrRep(irrep) ; i++ , moi++ )
			if ( !internal[n+i] )
			{
				mapRealToContinuous[moi-1] = mo;
				mapContinuousToReal[mo-1] = moi;
				mo++;
			}

		n += moirreps.getInIrRep(irrep);
	}
}

MOMapping::MOMapping(const MOMapping &map)
{
	maxMO = map.getMaxMO();
	mapRealToContinuous = new MOType[maxMO];
	mapContinuousToReal = new MOType[maxMO];
	for ( INT i=0 ; i<maxMO ; i++ )
	{
		mapRealToContinuous[i] = map.getContinuous(i+1);
		mapContinuousToReal[i] = map.getReal(i+1);
	}
}

MOMapping &	MOMapping::operator = (const MOMapping &map)
{
	if ( mapRealToContinuous )
		delete mapRealToContinuous;
	if ( mapContinuousToReal )
		delete mapContinuousToReal;
	maxMO = map.getMaxMO();
	mapRealToContinuous = new MOType[maxMO];
	mapContinuousToReal = new MOType[maxMO];
	for ( INT i=0 ; i<maxMO ; i++ )
	{
		mapRealToContinuous[i] = map.getContinuous(i+1);
		mapContinuousToReal[i] = map.getReal(i+1);
	}
	return *this;
}

MOMapping::~MOMapping()
{
	if ( mapRealToContinuous )
		delete[] mapRealToContinuous;
	if ( mapContinuousToReal )
		delete[] mapContinuousToReal;
}


void	MOMapping::setInternal(MRMOs *mrmos) const
{
	for ( INT i=0 ; i<maxMO ; i++ )
		if ( mapRealToContinuous[i]!=i+1 )
			mrmos->setInternal(i+1);
}

template <class Type>
void	MOMapping::reorderNoHoles(
	NExternalsDiag * &confTree, CIVectors<Type> * &ev,
	const MRMOs *mrmos) const
{


NExternalsSet	*confTreeNew = new NExternalsSet(
		mrmos,
		confTree->getNumberOfElectrons(),
		confTree->getMultiplicity(),
		confTree->getNumberOfRoots(),
		confTree->getTotalSymmetry());
			
//	cout << *mrmos << endl;
	
AVLMap<Configuration<MOType>, INT>	map(0);
	for ( MRTreeIterator i = confTree->firstInTree() ;
			!confTree->isLastInTree(i) ; confTree->nextInTree(i))
	{
	ConfigurationSAFNr<MOType>	conf(confTree->getConfigurationSAFNr(i));
//		cout << conf.getSAFNr() << " " << conf << endl;
		map[conf] = conf.getSAFNr();
		confTreeNew->add(conf);
	}
	
/*	cout << "===============================================" << endl;
	cout << "===============================================" << endl;
	cout << "===============================================" << endl;
	cout << "===============================================" << endl;

	for ( MRTreeIterator i = confTreeNew->firstInTree() ;
			!confTreeNew->isLastInTree(i) ; confTreeNew->nextInTree(i))
	{
		cout << confTreeNew->getConfigurationSAFNr(i) << endl;
	}
*/

	delete confTree;
	confTree = new NExternalsDiag(*confTreeNew);

CIVectors<Type>	*evNew = new CIVectors<Type>(
		ev->getDim(), ev->getN());

/*	cout << "===============================================" << endl;
	cout << "===============================================" << endl;
	cout << "===============================================" << endl;
	cout << "===============================================" << endl;
*/
	for ( MRTreeIterator i = confTree->firstInTree() ;
			!confTree->isLastInTree(i) ; confTree->nextInTree(i))
	{
	ConfigurationSAFNr<MOType>	conf(confTree->getConfigurationSAFNr(i));
	
	INT	safs = confTree->getNumberOfSpinAdaptedFunctions(
		conf.getNumberOfOpenShells());
//		cout << j2 << " " << j1 << (*ev)(j1) << " " << conf << endl;
		for ( INT k=0 ; k<ev->getN() ; k++ )
		{
		INT	j1 = map[conf];
		INT	j2 = conf.getSAFNr();
			for ( INT j=0 ; j<safs ; j++ )
				(*evNew)(j2++, k) = (*ev)(j1++, k);
		}
	}

	delete ev;
	ev = evNew;
	delete confTreeNew;
	
}



void	MOMapping::writeToStream(ostream &s)
{
	s << maxMO << endl;
	for ( INT i=0 ; i<maxMO ; i++ )
		s << mapRealToContinuous[i] << " ";
	s << endl;
	for ( INT i=0 ; i<maxMO ; i++ )
		s << mapContinuousToReal[i] << " ";
	s << endl;
}


ostream& operator<<(ostream & s, const MOMapping & map)
{
	if ( map.maxMO )
	{
		s << "MO mapping:" << endl;
		s << "real      : ";
		for ( INT i=0 ; i<map.maxMO ; i++ )
			s << setw(4) << i+1;
		s << endl;
		s << "continuous: ";
		for ( INT i=0 ; i<map.maxMO ; i++ )
			s << setw(4) << map.getContinuous(i+1);
		s << endl;
		s << endl;
		s << "continuous: ";
		for ( INT i=0 ; i<map.maxMO ; i++ )
			s << setw(4) << i+1;
		s << endl;
		s << "real      : ";
		for ( INT i=0 ; i<map.maxMO ; i++ )
			s << setw(4) << map.getReal(i+1);
		s << endl;
		s << endl;
	}
	else
		s << "no MO mapping." << endl;
	return s;
}

template void	MOMapping::reorderNoHoles(NExternalsDiag * &, CIVectors<double> * &,
		const MRMOs *) const;
template void	MOMapping::reorderNoHoles(NExternalsDiag * &, CIVectors<float> * &,
		const MRMOs *) const;

#include "../../Container/AVLMap.cc"

template class AVLMap<Configuration<MOType>, INT>;
template <> Configuration<MOType>*   AVLMap<Configuration<MOType>, INT>::_target_item = NULL;     // add/del_item target
template <> AVLNode<Configuration<MOType>, INT>* 
	AVLMap<Configuration<MOType>, INT>::_found_node = NULL; // returned added/deleted node

#include "../../Container/Map.cc"
template class Map<Configuration<MOType>, INT>;

