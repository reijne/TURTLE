//***********************************************************************
//
//	Name:			EnergyMap.h
//
//	Description:	
//
//	Author:			Michael Hanrath
//
//	Version:		0.0
//
//	Date:			19.04.1997
//
//
//
//
//***********************************************************************


#ifndef __EnergyMap_h
#define __EnergyMap_h

#include "../../../config.h"


#include "../../../lib/Container/AVLMap.h"

#include "../../../lib/QM/Configuration/Configuration.h"
#include "../../../lib/QM/MRTree/EnergyMap/EnergyEntry.h"
#include "../../../lib/QM/MO/Iterators/MOList.h"

class TupelStructureSel;
class MORestriction;
class MOStatistics;
class NExternalsSet;
class PTSum;
template <class T> class Histogram;
template <class RepMatType, class VectorType> class EnlargeReferenceSpace;

class EnergyMap {
public:
enum	EstimationMode	{ EpsteinNesbet, Wenzel };


typedef double	RepMatType;
typedef double	VectorType;

	EnergyMap(INT roots, EstimationMode estimationMode);
	~EnergyMap();


	
	void	test();
	
	INT	getNumberOfConfigurations() const;
	
	INT	getEntries(
		const Configuration<MOType> &,
		INT safs,
		EnergyType * & nominator,
		EnergyType * & denominator
		);

	void	select(
		TupelStructureSel *tupel,
		const MORestriction *moRestriction,
		const MOStatistics *moStatistics,
		EnlargeReferenceSpace<RepMatType, VectorType> *enlargeRef,
		EnergyType threshold,
		INT	randThresh,
		Histogram<EnergyType> **hist,
		PTSum *ptTotalSum,
		PTSum **ptSum,
		INT	&Confs,
		INT &CSFs,
		INT	&RestrictedConfs,
		INT &RestrictedCSFs,
		EnergyType *RestrictedESum,
		INT	selectAll,
		EnergyEntry::SelectionMode mode = EnergyEntry::AmountSquareOfSum
		);
	
	void	select(
		TupelStructureSel *tupel,
		NExternalsSet *preSelConfs, INT nExt,
		EnergyType *eSum,
		EnergyEntry::SelectionMode mode = EnergyEntry::AmountSquareOfSum
		);
	
	
	
private:
INT	nRoots;					// number of roots


AVLMap<Configuration<MOType>, EnergyEntry *>
	*keyAVLMap;				// AVLMap implementation of GNU C++ Library
	
EstimationMode	estimationMode;
};



inline
EnergyMap::EnergyMap(INT roots, EstimationMode _estimationMode)
{
	keyAVLMap = new AVLMap<Configuration<MOType>, EnergyEntry *>(NULL);
	nRoots = roots;
	estimationMode = _estimationMode;
}

inline
EnergyMap::~EnergyMap()
{
Pix	i = keyAVLMap->first();
	while ( i )
	{
		delete keyAVLMap->contents(i);
		keyAVLMap->next(i);
	}
	delete keyAVLMap;
}


inline
INT	EnergyMap::getNumberOfConfigurations() const
{	return keyAVLMap->length();	}

inline
INT	EnergyMap::getEntries(
	const Configuration<MOType> & key,
	INT safs,
	EnergyType * & nominator,
	EnergyType * & denominator
)
{
EnergyEntry	**p = &(*keyAVLMap)[key];
	if ( *p )
	{
		nominator = (*p)->getNominatorP();
		denominator = (*p)->getDenominatorP();
		return 1;
	}
	
	*p = new EnergyEntry(nRoots, safs);
	nominator = (*p)->getNominatorP();
	denominator = (*p)->getDenominatorP();
	return 0;
}

inline
void	EnergyMap::test()
{	cout << "EnergyMap::roots=" << nRoots << endl;	}

#endif
