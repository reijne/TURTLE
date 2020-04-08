//***********************************************************************
//
//	Name:			EnergyMap.cc
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




#include "EnergyMap.h"

#include "../../../lib/QM/MRTree/Set/NExternalsSet.h"
#include "../../../lib/QM/MRTree/Set/InternalConfsSet.h"
#include "../../../lib/QM/MRTree/Set/TupelStructureSet.h"
#include "../../../lib/QM/MRTree/Set/extMOsSet.h"
#include "../../../lib/QM/MRTree/Sel/TupelStructureSel.h"
#include "../../../lib/QM/MRTree/Sel/extMOsSel.h"
#include "../../../lib/QM/MRTree/Base/extEntry.h"
#include "../../../lib/Math/etc/Histogram.h"
#include "../../../lib/QM/MRTree/EnergyMap/PTSum.h"
#include "../../../lib/QM/MO/Iterators/MORestrictionState.h"
#include "../../../lib/QM/MO/MOStatistics.h"

#include "EnlargeReferenceSpace.h"




void	EnergyMap::select(
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
	EnergyEntry::SelectionMode mode)
{
Pix	i = keyAVLMap->first();	
EnergyEntry	*entry;

//	cout << moRestriction << endl;
//	cout << moRestriction->getRestriction() << endl;
MORestrictionState	*restrictionState = NULL;

	if ( moRestriction )
		restrictionState = new MORestrictionState(*moRestriction, *tupel);

//	cout << "AAAAAAAAA" << endl;

//	cout << *tupel << endl;

	while ( i )
	{
		entry = keyAVLMap->contents(i);
	Configuration<MOType>	*p = &keyAVLMap->key(i);
	INT	restricted = 0;

		if ( moRestriction )
		{
			MORestrictionState	restrictionState2(*moRestriction, *p);
			restrictionState2 |= *restrictionState;
			if ( !restrictionState2.check() )
				restricted = 1;
		}


/*	Configuration<MOType> conf(*tupel);
		conf += *p;
		if ( moRestriction->check(conf)!=restrictionState2.check() )
			cout << restrictionState2.check() << " " << conf << endl;
*/
//		cout << *p << endl;
	INT	SAFs;
	extMOsSel	*pSel = tupel->operator [] (p->getNumberOfOpenShells()/2);

		SAFs = pSel->getSAFInc();
		if ( !restricted )
		{
			Confs++;
			CSFs += SAFs;


			switch ( estimationMode ) {
			case EpsteinNesbet:
				{
				EnergyEntry::Status	flag = entry->select(threshold, mode);

//-----------------------------------------------------------------
//	random selection
					if ( !(flag || selectAll) && randThresh )
						if ( rand()>randThresh )
							flag = EnergyEntry::Random;
//-----------------------------------------------------------------


					if ( flag || selectAll )
					{
						pSel->add(
							*p,
							nRoots, entry->getEnergyP(),
							SAFs, entry->getCICoefP(),
								(flag == EnergyEntry::SelectionFlagserated && !selectAll)
								| ((selectAll>0) << 1)
							);
			//				cout << "OOOOOOOOOOOO " << keyAVLMap->key(i) << endl;
					}

					ptTotalSum->count(entry->getEnergyP(), SAFs, 
						selectAll || (flag==EnergyEntry::Random));

					for ( INT j=0 ; j<nRoots ; j++ )
					{
						hist[j]->count(entry->getEnergy(j));
						ptSum[j]->count(entry->getEnergy(j), SAFs, 
							selectAll || (flag==EnergyEntry::Random));
					}
				}
				break;


			case Wenzel:
				{
				Configuration<MOType> conf(*tupel);
				const MOType *pmo = p->getOpenShellP();
					for ( INT j=0 ; j<p->getNumberOfOpenShells() ; j++ )
						conf.insertOpenMO(pmo[j]);
					pmo = p->getClosedShellP();
					for ( INT j=0 ; j<p->getNumberOfClosedShells() ; j++ )
						conf.insertClosedMO(pmo[j]);
				const EnergyType *pE = enlargeRef->calcEnlargementEnergy(conf);


				INT	flag = (fabs(pE[0])>threshold);

					if ( flag || selectAll )
					{
						pSel->add(
							*p,
							nRoots, pE,
							SAFs, entry->getCICoefP(),
								(flag == EnergyEntry::SelectionFlagserated && !selectAll)
								| ((selectAll>0) << 1)
							);
					}

					ptTotalSum->count(pE, SAFs, selectAll);

					for ( INT j=0 ; j<nRoots ; j++ )
					{
						hist[j]->count(pE[j]);
						ptSum[j]->count(pE[j], SAFs, selectAll);
					}
				}
				break;
			}
		}
		else
		{
			RestrictedConfs++;
			RestrictedCSFs += SAFs;
		
			switch ( estimationMode ) {
			case EpsteinNesbet:
				entry->select(threshold, mode);
				for ( INT j=0 ; j<nRoots ; j++ )
					RestrictedESum[j] += entry->getEnergyP()[j];
				break;


			case Wenzel:
				{
				Configuration<MOType> conf(*tupel);
				const MOType *pmo = p->getOpenShellP();
					for ( INT j=0 ; j<p->getNumberOfOpenShells() ; j++ )
						conf.insertOpenMO(pmo[j]);
					pmo = p->getClosedShellP();
					for ( INT j=0 ; j<p->getNumberOfClosedShells() ; j++ )
						conf.insertClosedMO(pmo[j]);

				const EnergyType *pE = enlargeRef->calcEnlargementEnergy(conf);
					for ( INT j=0 ; j<nRoots ; j++ )
						RestrictedESum[j] += pE[j];
				}
				break;
			}
		}
		keyAVLMap->next(i);
	}
	if ( restrictionState )
		delete restrictionState;
//	cout << "BBBBBBBBBB" << endl;
}




void	EnergyMap::select(
	TupelStructureSel *tupel,
	NExternalsSet *preSelConfs, INT nExt,
	EnergyType	*eSum,
	EnergyEntry::SelectionMode mode)
{
	InternalConfsSet	*t = new InternalConfsSet(nExt, NULL);
Pix	iExt = preSelConfs->seek(t);
	delete t;
TupelStructureSet	*preSelTupel = NULL;
Pix	ind = NULL;
	if ( iExt )
	{
	InternalConfsSet	*intConfs = (*preSelConfs)[iExt];
	TupelStructureSet	*indTupel = new TupelStructureSet(NULL);
		*((TupelStructureBase *) indTupel) = *((TupelStructureBase *) tupel);
//	cout << *indTupel << endl;
	
		ind = intConfs->seek(indTupel);
		delete indTupel;

		if ( ind ) 
			preSelTupel = (*intConfs)[ind];
	}
//	cout << "ind=" << ind << endl;

Pix	i = keyAVLMap->first();	
EnergyEntry	*entry;

	while ( i )
	{
		entry = keyAVLMap->contents(i);
	INT	sumUp = 1;
		if ( ind && entry )
		{
		Configuration<MOType>	*p = &keyAVLMap->key(i);
		extMOsSet	*t = new extMOsSet(
			p->getNumberOfOpenShells(), p->getNumberOfClosedShells());
		Pix pp = preSelTupel->seek(t);
			if ( pp )
			{
			extMOsSet	*ppreSel = preSelTupel->operator [] (pp);
				if ( ppreSel )
				{
				extEntry	*ext = new extEntry(ppreSel, 0, NULL, 0, NULL, *p, 0);
	//				cout << "---------------------" << endl;
	//				ppreSel->writeToStream(cout);
	//				cout << endl;
	//				cout << "---------------------" << endl;
					sumUp = !ppreSel->contains(ext);
	//				cout << "ext=" << *ext << endl;
	//				cout << "sumUp=" << sumUp << endl;
					delete ext;
				}
			}
			delete t;
		}
		if ( sumUp )
		{
			entry->select(0.0, mode);
			for ( INT j=0 ; j<nRoots ; j++ )
				eSum[j] += entry->getEnergy(j);
		}
		keyAVLMap->next(i);
	}
}
