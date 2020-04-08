//***********************************************************************
//
//	Name:			MatchingMOListIterator.cc
//
//	Description:	iterators for MO lists
//
//	Author:			Michael Hanrath
//
//	Version:		0.0
//
//	Date:			26.03.1997
//
//
//
//
//
//***********************************************************************

#include "MatchingMOListIterator.h"

#include "TupelIterator.h"
#include "OrderMatrixIterator.h"


#include <stdlib.h>


MatchingMOListIterator::MatchingMOListIterator(
	ConfigurationSAFNr<MOType> base,
	Configuration<MOType> toMatch,
	const MOList *_molist,
	INT minOpenShells)
{
	molist = _molist;
	matchList = NULL;
	SLIter = NULL;
const MRMOs	*mrmos = molist->getMRMOs();

Configuration<MOType>	extBase = base-base.getInternal();

//	cout << "base: " << base << endl;
//	cout << "base (internal): " << base.getInternal() << endl;
//	cout << "extBase: " << extBase << endl;
//	cout << "toMatch: " << toMatch << endl;
//	cout << "molist: " << *molist << endl;

INT	NintDif = base.getInternal().getNumberOfElectrons() - 
		toMatch.getNumberOfElectrons();
			

//	cout << "NintDif=" << NintDif << endl;


	if ( abs(NintDif)>2 )
	{
//		cout << "EXIT." << endl;
		return;
	}
		


	if ( Configuration<MOType>::calcExcitationOrder(base.getInternal(), toMatch)>2 )
	{
//		cout << "EXIT." << endl;
		return;
	}
		
INT	dif = Configuration<MOType>::calcExcitationOrder(base, toMatch);// + molist->getNumberOfExternals();
	if ( NintDif>0 )
		dif -= NintDif;
//	cout << "dif=" << dif << endl;
	
	
	
//	cout << "molist->getNumberOfExternals(): " << molist->getNumberOfExternals() << endl;


	//	diminish total excitation by possible external matches
	
Configuration<MOType>	extBaseMatch = extBase;	
	for ( INT i=0 ; i<molist->getNumberOfExternals() ; i++ )
	{
	IrRep ir = molist->getIrRepExt(i);
	INT	flag = 0;
//		cout << "ir=" << ir << ", extBaseMatch=" << extBaseMatch << ", dif=" << dif << endl;
		for ( INT j=0 ; j<extBaseMatch.getNumberOfOpenShells() ; j++ )
			if ( ir==mrmos->getIrRep(extBaseMatch.getOpenShell(j)) )
			{
				dif--;
				extBaseMatch.annihilate(extBaseMatch.getOpenShell(j));
				flag = 1;
				break;
			}
		if ( flag )
			continue;

		for ( INT j=0 ; j<extBaseMatch.getNumberOfClosedShells() ; j++ )
			if ( ir==mrmos->getIrRep(extBaseMatch.getClosedShell(j)) )
			{
				dif--;
				extBaseMatch.annihilate(extBaseMatch.getClosedShell(j));
				break;
			}
	}
	
//	cout << "dif=" << dif << endl;
	if ( dif>2 )
	{
//		cout << "EXIT." << endl;
		return;
	}

//	cout << "START" << endl;


//------------------------------------------------------------------------

INT	 extIrReps = molist->getNumberOfExternalIrReps();
INT* nBaseMOsInIrRep = new INT[extIrReps];
INT* nBaseElectronsInIrRep = new INT[extIrReps];
MOType* *baseMOsInIrRep = new MOType*[extIrReps];

INT* iDif = new INT[extIrReps+1];

	for ( INT i=0 ; i<extIrReps ; i++ )
	{
	INT	k = 0;
		nBaseMOsInIrRep[i] = 0;
		nBaseElectronsInIrRep[i] = 0;
		baseMOsInIrRep[i] = new MOType[base.getNExt()];
	IrRep ir = molist->getExternalIrRep(i);

		for ( INT j=0 ; j<extBase.getNumberOfOpenShells() ; j++ )
			if ( ir==mrmos->getIrRep(extBase.getOpenShell(j)) )
			{
				nBaseMOsInIrRep[i]++;
				nBaseElectronsInIrRep[i] += 1;
				baseMOsInIrRep[i][k++] = extBase.getOpenShell(j);
			}

		for ( INT j=0 ; j<extBase.getNumberOfClosedShells() ; j++ )
			if ( ir==mrmos->getIrRep(extBase.getClosedShell(j)) )
			{
				nBaseMOsInIrRep[i]++;
				nBaseElectronsInIrRep[i] += 2;
				baseMOsInIrRep[i][k++] = extBase.getClosedShell(j);
			}
	}

	iDif[0] = dif;

//	if ( NintDif<0 )
//		iDif[0] += NintDif;

//	cout << "iDif[0]=" << dif << endl;
	
//	cout << "extIrReps=" << extIrReps << endl;

	dcGenBase.calcDiffConf(base, toMatch);

	matchList = new SLList<MatchingMOList>;
	if ( extIrReps>0 )
	{
	//	prepare tupel structure iterators
TupelIterator* *tupels = new TupelIterator*[extIrReps];
	for ( INT i=0 ; i<extIrReps ; i++ )
	{
		tupels[i] = new TupelIterator(molist->getInExternalIrRep(i));
//		cout << "inextirrep=" << molist->getInExternalIrRep(i) << endl;
	}


	//	****************** loops over tupel structures
	{
INT	i = 0;
	tupels[i]->first();
INT	nOpen = 0;

	for ( ; ; )
	{
		nOpen += tupels[i]->getNumberOfOpenShells();
//		cout << *tupels[i] << endl;
		if ( i<extIrReps-1 )
		{
			tupels[++i]->first();
			continue;
		}
//		cout << "nOpen= " << nOpen << ", minOpenShells= " << minOpenShells << endl;
		if ( i==extIrReps-1 && nOpen>=minOpenShells )
		{
			//	prepare order matrix iterators
		OrderMatrixIterator* *orders = new OrderMatrixIterator*[extIrReps];
			for ( INT j=0 ; j<extIrReps ; j++ )
				orders[j] = new OrderMatrixIterator(
					nBaseMOsInIrRep[j], tupels[j]->getNumberOfShells());
				
				
			//	################## loops over order matrices
			{
		INT	j = 0;
			orders[j]->first();
			for ( ; ; )
			{
//				cout << *orders[j] << endl;
				iDif[j+1] = iDif[j];
//				cout << "v. iDif[" << j+1 << "]=" << iDif[j+1] << endl;
				{
				INT	k = 0;
					iDif[j+1] += abs(nBaseElectronsInIrRep[j] -
						(tupels[j]->getNumberOfOpenShells() + 
						2*tupels[j]->getNumberOfClosedShells()));
					if ( nBaseMOsInIrRep[j]>0 )
					{
					INT	d = nBaseElectronsInIrRep[j] - nBaseMOsInIrRep[j];
						for ( ; k<tupels[j]->getNumberOfOpenShells() ; k++ )
							iDif[j+1] += orders[j]->isEqualInColumn(k) ? 0 : 1;
						for ( ; k<tupels[j]->getNumberOfShells() ; k++ )
							iDif[j+1] += orders[j]->isEqualInColumn(k) ? ((d--) ? 0 : 1) : 2;
					}
				}
//				cout << "n. iDif[" << j+1 << "]=" << iDif[j+1] << endl;
				if ( iDif[j+1] <= 2 )
				{
					if ( j<extIrReps-1 )
					{
						orders[++j]->first();
						continue;
					}

					if ( j==extIrReps-1 )
					{
					MatchingMOList	match;
						match.moiter = new MOListIterator;
						match.dcGen = new DiffConf<GeneralizedMO>(dcGenBase);
						match.extNotRunning = new Configuration<MOType>;


						//	the following works for max. 2 external MOs
						//	on right hand side
						//	to generalize to more than 2 external MOs 
						//	on right hand side a selection of 2 
						//	out of n external MOs has to be performed

					INT	sig = 0;
					INT	empty = 0;
						for ( INT k=0 ; k<=extIrReps-1 ; k++ )
						{
						IrRep	irrep = molist->getExternalIrRep(k);
//							cout << "####:" << 
//								mrmos->getIrRepIntExtStart(
//									irrep, 1) << ", " <<
//								mrmos->getIrRepIntExtEnd(
//									irrep, 1) << endl;

							empty = orders[k]->calcMOListIterator(
								mrmos->getIrRepIntExtStart(
									irrep, 1),
								mrmos->getIrRepIntExtEnd(
									irrep, 1),
								tupels[j]->getNumberOfClosedShells(),
								baseMOsInIrRep[k],
								*match.moiter,
								*match.extNotRunning,
								*match.dcGen,
								sig,
								irrep);

							if ( empty )
								break;

//							cout << "### moiter= " << *match.moiter << endl;
//							cout << "### dcGen= " << *match.dcGen << endl;
						}


//						cout << "moiter= " << *match.moiter << endl;
						
//						cout << "dcGen= " << *match.dcGen << endl;
//						cout << "extNotRunning= " << *match.extNotRunning << endl;
						*match.extNotRunning += match.dcGen->simplify();
//						cout << "extNotRunning simp.= " << *match.extNotRunning << endl;
//						cout << "dcGen simp.= " << *match.dcGen << endl;


						if ( empty )
						{
							delete match.moiter;
							delete match.dcGen;
							delete match.extNotRunning;
						}
						else
							matchList->append(match);
					}
				}
				orders[j]->next();
				while ( j>=0 && orders[j]->isEnd() )
				{
					if ( --j<0 )
						break;
					orders[j]->next();
				}

				if ( j<0 )
					break;
			}
			}
			//	################## end loops over order matrices


			//	release order matrix iterators
			for ( INT j=0 ; j<extIrReps ; j++ )
				delete orders[j];
			delete[] orders;


		}
	
//		cout << "i=" << i << endl;

		nOpen -= tupels[i]->getNumberOfOpenShells();
		tupels[i]->next();
		while ( i>=0 && tupels[i]->isEnd() )
		{
			if ( --i<0 )
				break;
			nOpen -= tupels[i]->getNumberOfOpenShells();
			tupels[i]->next();
//		cout << "i=" << i << endl;
		}

		if ( i<0 )
			break;
	}
	}
	//	****************** end loops over tupel structures




	//	release tupel structure iterators
	for ( INT i=0 ; i<extIrReps ; i++ )
		delete tupels[i];


	for ( INT i=0 ; i<extIrReps ; i++ )
		delete[] baseMOsInIrRep[i];

	delete[] tupels;

	}	// if ( extIrReps>0 )
	else
	{
	MatchingMOList	match;
		match.moiter = new MOListIterator(0, 0, 0);
		match.dcGen = new DiffConf<GeneralizedMO>(dcGenBase);
		match.extNotRunning = new Configuration<MOType>;

		matchList->append(match);
	}

        delete[] iDif;
        delete[] baseMOsInIrRep;
        delete[] nBaseElectronsInIrRep;
        delete[] nBaseMOsInIrRep;

	SLIter = matchList->first();
//	printf("SLIter = %lx\n", SLIter);

//	printf("%lx %lx\n",(* matchList)(SLIter).dcGen, (*matchList)(SLIter).moiter);
//	cout << *(*matchList)(SLIter).dcGen << endl;

	
//	cout << "END" << endl;
}


MatchingMOListIterator::~MatchingMOListIterator()
{
	if ( matchList )
	{
		SLIter = matchList->first();
		while ( SLIter )
		{
			delete (*matchList)(SLIter).moiter;
			delete (*matchList)(SLIter).dcGen;
			delete (*matchList)(SLIter).extNotRunning;
			matchList->next(SLIter);
		}
		delete matchList;
	}
}


const DiffConf<GeneralizedMO> &	MatchingMOListIterator::getGenDiffConf() const
{	return *(*matchList)(SLIter).dcGen;	}

const MOListIterator &	MatchingMOListIterator::getMOListIterator() const
{	return *(*matchList)(SLIter).moiter;	}

const Configuration<MOType> &	MatchingMOListIterator::getExtNotRunning() const
{	return *(*matchList)(SLIter).extNotRunning;	}

void	MatchingMOListIterator::next()
{	matchList->next(SLIter);	}

INT	MatchingMOListIterator::isEnd() const
{	return SLIter==0;	}


#include "../../../Container/SLList.cc"
template class SLList<MatchingMOList>;

