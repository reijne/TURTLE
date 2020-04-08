//***********************************************************************
//
//	Name:			Selector.cc
//
//	Description:	
//
//	Author:			Michael Hanrath
//
//	Version:		0.0
//
//	Date:			23.03.1997
//
//
//
//
//
//***********************************************************************

using namespace std;

#include "Selector.h"
#include "../../../lib/QM/IO/Fortran/Fort31FirstRecord.h"
#include "../../../lib/QM/MRTree/Sel/InternalConfsSel.h"
#include "../../../lib/QM/MO/Iterators/MOList.h"
#include "../../../lib/QM/MO/Iterators/MOIterator.h"


#include "../../../lib/QM/IntegralIndex/TwoElectronIntegralTriadeIndex.h"
#include "../../../lib/QM/IntegralIndex/TwoElectronIntegralIndexTemplate.h"

#include "../../../lib/QM/Cache/TableKey.h"
#include "../../../lib/QM/Cache/VarSizeReadOnlyCache.h"
#include "EnlargeReferenceSpace.h"

#include "../../../lib/QM/MO/Iterators/MatchingMOListIterator.h"

#include "../../../lib/QM/RepresentationMatrices/PTReferenceSpace.h"
#include "../../../lib/QM/RepresentationMatrices/RepDiag.h"
#include "../../../lib/QM/IO/TimeTicks.h"

#include "EnergyMap.h"
#include "../../../lib/QM/MRTree/EnergyMap/PTSum.h"

#include "../../../lib/Math/etc/Histogram.h"

#include "../../../lib/QM/Parallel/SharedMemory.h"

#include "../../../lib/QM/IO/Verbosity.h"

#include <iomanip>

template <class MatrixType, class VectorType>
Selector<MatrixType, VectorType>::Selector(MRConfInput _mrinp, NExternalsSet *_preSelConfs,
		INT cacheEntries, INT cacheMemory):
		CICalculation<MatrixType, VectorType>(_mrinp.getMultiplicity(),
		0, *_mrinp.getMRMOs(), _mrinp.getMOIntegralFile(), this->storeAll)
{
	mrinp = _mrinp;
	tablecaseGen = new TableCase<GeneralizedMO>(this->binom);
			
	cache = new VarSizeReadOnlyCache<TableKey, HMatElements<MatrixType, VectorType> >
			(cacheEntries, cacheMemory);
	
	preSelConfs = _preSelConfs;
	if ( preSelConfs )
	{
		preESum = new EnergyType[mrinp.getNumberOfRoots()];
		memset(preESum, 0, mrinp.getNumberOfRoots()*sizeof(EnergyType));
	}
	else
		preESum = NULL;

	GeneralizedMO::mrmos = &this->mrmos;



	minOpenShells = mrinp.getMultiplicity()-1;



//	cout << "calculating P=5 representation matrices..." << endl;
	repMatsP5 = new HMatElements<MatrixType, VectorType> * [MAXSGATABLEOPENSHELLS+1];
	memset(repMatsP5, 0, (MAXSGATABLEOPENSHELLS+1)*sizeof(HMatElements<MatrixType, VectorType> *));
	
	for ( INT i=minOpenShells ; i<=MAXSGATABLEOPENSHELLS ; i+=2 )
	{
//		cout << i << " open shells:" << endl;
		this->tablecase->setNr(i, 28);	// P=5, case 14
		repMatsP5[i] =  new HMatElements <MatrixType, VectorType>
			(*(*cache)[TableKey(*this->tablecase)]);
//		cout << *repMatsP5[i] << endl;
	}

//==============================================================================
//==============================================================================
//==============================================================================

	genConfs = new NExternalsSel(
		mrinp.getMRMOs(),
		mrinp.getNumberOfElectrons(),
		mrinp.getMultiplicity(),
		mrinp.getExcitationLevel(),
		mrinp.getNumberOfRoots(),
		mrinp.getRefConfSet()
		);

	genConfCount = new INT[mrinp.getExcitationLevel()+1];
	memset(genConfCount, 0, (mrinp.getExcitationLevel()+1)*sizeof(INT));
	genCSFCount = new INT[mrinp.getExcitationLevel()+1];
	memset(genCSFCount, 0, (mrinp.getExcitationLevel()+1)*sizeof(INT));

	genRestrictedConfCount = new INT[mrinp.getExcitationLevel()+1];
	memset(genRestrictedConfCount, 0, (mrinp.getExcitationLevel()+1)*sizeof(INT));
	genRestrictedCSFCount = new INT[mrinp.getExcitationLevel()+1];
	memset(genRestrictedCSFCount, 0, (mrinp.getExcitationLevel()+1)*sizeof(INT));

	RestrictedESum = new EnergyType[mrinp.getNumberOfRoots()];
	memset(RestrictedESum, 0, (mrinp.getNumberOfRoots())*sizeof(EnergyType));

	enlargeRef = NULL;
	if ( mrinp.getEstimationMode()==EnergyMap::Wenzel )
		enlargeRef = new EnlargeReferenceSpace<MatrixType, VectorType>(mrinp);
	
	diagHist = new Histogram<EnergyType>(
//		enlargeRef->getRootEnergy(0), enlargeRef->getRootEnergy(0)+20,
		-100, -80,
		Histogram<EnergyType>::Linear, 100);
}



template <class MatrixType, class VectorType>
Selector<MatrixType, VectorType>::~Selector()
{
	delete diagHist;
	if ( enlargeRef )
		delete enlargeRef;

	delete genConfCount;
	delete genCSFCount;
	delete genRestrictedConfCount;
	delete genRestrictedCSFCount;
	delete RestrictedESum;

	delete genConfs;

	delete cache;

	delete tablecaseGen;
	for ( INT i=0 ; i<MAXSGATABLEOPENSHELLS ; i++ )
		if ( repMatsP5[i] )
			delete repMatsP5[i];
	delete repMatsP5;
	if ( preSelConfs )
		delete preSelConfs;
	if ( preESum )
		delete preESum;
}




template <class MatrixType, class VectorType>
void	Selector<MatrixType, VectorType>::useAnnihilators()
{
TimeTicks	ticks;

	cout << endl << "using annihilators... " << flush;
	ticks.start();

	for ( INT i=1 ; i<genConfs->getNumberOfElements() ; i++ )
		(*genConfs)[i] = new InternalConfsSel (i, genConfs);

	
	
	// use internal annihilators
ConfigurationSet	&references = mrinp.getRefConfSet();
Pix	iRef = references.first();
	while ( iRef )
	{
		for ( INT i=0 ; i<genConfs->getNumberOfElements() ; i++ )
		{
		MOIterator	internalAnnihilators(i, &this->mrmos, 0);
			while ( !internalAnnihilators.isEnd() )
			{
			Configuration<MOType>	conf(references(iRef));
				conf -= internalAnnihilators;
				if ( conf.getNumberOfElectrons() )
				{
				TupelStructureSel *tupel = new TupelStructureSel(
						i, conf.calcIrRep(this->mrmos),
						conf, (*genConfs)[i]);

					if ( !(*genConfs)[i]->contains(tupel) )
						(*genConfs)[i]->add(tupel);
					else
						delete tupel;

				}
				internalAnnihilators.next();
			}
		}
		references.next(iRef);
	}
	ticks.stop();

	cout << ticks << endl;
}




template <class MatrixType, class VectorType>
void	Selector<MatrixType, VectorType>::useCreators(
	const PTReferenceSpace<MatrixType, VectorType> *_PTReference,
	EnergyType Threshold)
{
TimeTicks	ticks;
	cout << endl << "using creators... " << flush;
	ticks.start();

	PTReference = _PTReference;
	//	loop over number of creations
	for ( INT ncreate=1 ; ncreate<genConfs->getNumberOfElements() ; ncreate++ )
	{
InternalConfsSel	*intern = (*genConfs)[ncreate];

		//	loop over intern-i configurations
		for ( ContainerIterator iter = intern->first() ;
			!intern->isLast(iter) ; intern->next(iter) )
		{
	TupelStructureSel *tupel = intern->operator [] (iter);


//		cout << "====================" << endl;
//		cout << *((Configuration<MOType> *) tupel) << endl;

		IrRep	irrep[ncreate+1];
			irrep[0] = this->mrmos.getNumberOfIrReps() - 1;
			irrep[1] = 0;

		
//**************************************************************************
		INT	l = 1;
			//	loops over symmetries
			for ( ; ; irrep[l]++ )
			{
				if ( irrep[l] >= irrep[l-1]+1 )
				{
					l--;
					if ( l<1 )
						break;
					else
						continue;
				}
				if ( l<ncreate )
				{
					l++;
					irrep[l] = -1;
					continue;
				}
				else
				{
					//	check matching symmetry of total wave function
				IrRep prod = tupel->getIrRep();
					for ( INT i=0 ; i<ncreate ; i++ )
						prod = this->mrmos.getProd(prod, irrep[i+1]);

					if ( prod != genConfs->getTotalSymmetry() )
						continue;
						
				
					//	generate all possible internal/external combinations
				INT	isInternal[ncreate];
					for ( INT i=0 ; i<(2 << (ncreate-1)) ; i++ )
					{
						for ( INT j=0 ; j<ncreate ; j++ )
							isInternal[j] = (i & (1 << j)) ? 0 : 1;

						//	check for ambiguous combinations
						if ( checkAmbiguity(ncreate, isInternal, irrep+1) )
						{
						MOList	molist(ncreate, &this->mrmos, isInternal, irrep+1);
						
							intern->addToGenerated(molist.getTotalNumber());
//							cout << molist << endl;
							select(tupel, molist, Threshold);
//							cout << "...................." << endl;
						}
					}
				}
			}	//	loops over symmetries
//**************************************************************************
		}	//	loop over intern-i configurations
	}	//	loop over number of creations
	ticks.stop();
	cout << ticks << endl;
	
INT	confSum = 0;
INT	CSFSum = 0;
INT	confSumRestricted = 0;
INT	CSFSumRestricted = 0;
	cout << endl << endl;
	if ( mrinp.getMORestriction() )
	{
		cout << "total generated configurations/CSFs:" << endl << endl;
		cout << "              restrictions fulfilled     restrictions not fulfilled" << endl;
		cout << "                confs.        CSFs           confs.        CSFs" << endl;
		for ( INT i=0 ; i<genConfs->getNumberOfElements() ; i++ )
		{
			cout << "intern-" << i << ": "
				<< setw(12) << genConfCount[i]
				<< setw(12) << genCSFCount[i] << "     "
				<< setw(12) << genRestrictedConfCount[i]
				<< setw(12) << genRestrictedCSFCount[i] << endl;
			confSum += genConfCount[i];
			CSFSum += genCSFCount[i];
			confSumRestricted += genRestrictedConfCount[i];
			CSFSumRestricted += genRestrictedCSFCount[i];
		}
		cout << "          -------------------------------------------------------" << endl;
		cout << "total   : " << setw(12) << confSum
			<< setw(12) << CSFSum << "     "
			<< setw(12) << confSumRestricted 
			<< setw(12) << CSFSumRestricted << endl;
		cout << endl << endl;
	}
	else
	{
		cout << "total generated configurations/CSFs:" << endl << endl;
		cout << "                confs.        CSFs" << endl;
		for ( INT i=0 ; i<genConfs->getNumberOfElements() ; i++ )
		{
			cout << "intern-" << i << ": "
				<< setw(12) << genConfCount[i]
				<< setw(12) << genCSFCount[i] << endl;
			confSum += genConfCount[i];
			CSFSum += genCSFCount[i];
		}
		cout << "          ------------------------" << endl;
		cout << "total   : " << setw(12) << confSum
			<< setw(12) << CSFSum << endl;
		cout << endl << endl;
	}
}


template <class MatrixType, class VectorType>
INT	Selector<MatrixType, VectorType>::checkAmbiguity(INT n, INT *isInternal, IrRep *irrep) const
{
	for ( INT i=0 ; i<n-1 ; i++ )
		if ( irrep[i]==irrep[i+1] && isInternal[i]<isInternal[i+1] )
			return 0;
	return 1;
}




template <class MatrixType, class VectorType>
void	Selector<MatrixType, VectorType>::select(
	TupelStructureSel *tupel,
	const MOList & molist,
	EnergyType Threshold)
{

Configuration<MOType>* conf = new Configuration<MOType>[molist.getNumberOfInternals()+1];
		conf[0] = *tupel;
//		cout << "*******" << conf[0] << endl;
//		cout << molist << endl;

INT	selectInternal = mrinp.getSelectInternal();

INT	selectNthExcitation = mrinp.getSelectNthExcitation();


TwoElectronIntegralIndexTemplate	CoulombIndexT, ExchangeIndexT;
TwoElectronIntegralIndex<GeneralizedMO>	CoulombIndexGen, ExchangeIndexGen;
	
INT	j[molist.getNumberOfInternals()+1];
	j[0] = this->mrmos.getIrRepIntExtStart(molist.getIrRep(0), 0);
INT	i = 0;
	//	create internal MOs
	for ( ; ; j[i]++ )
	{
		if ( molist.getNumberOfInternals()>0 )
		{
			if ( j[i] > this->mrmos.getIrRepIntExtEnd(molist.getIrRep(i), 0) )
			{
				i--;
				if ( i<0 )
					break;
				else
					continue;
			}
			conf[i+1] = conf[i];
			if ( conf[i+1].create(j[i]) )
				continue;
		}
		else
			i = -1;
		
//			cout << i << ": " << conf[i+1] << endl;


		if ( i<molist.getNumberOfInternals()-1 )
		{
			i++;
			j[i] = this->mrmos.getIrRepIntExtStart(molist.getIrRep(i), 0) - 1;
			continue;
		}
		else
		{
		TupelStructureSel *tupelWInt;
		
	
			if ( molist.getNumberOfInternals()>0 )
			{
				tupelWInt = new TupelStructureSel(
					molist.getNumberOfExternals(),
					conf[i+1].calcIrRep(this->mrmos),
					conf[i+1], (*genConfs)[molist.getNumberOfExternals()]);

				if ( !(*genConfs)[molist.getNumberOfExternals()]->contains(tupelWInt) )
					(*genConfs)[molist.getNumberOfExternals()]->add(tupelWInt);
				else
				{
					delete tupelWInt;
					goto ende;
				}
			}
			else
				tupelWInt = tupel;
			

			INT	maxOpenShells = conf[i+1].getNumberOfOpenShells() +
				molist.getNumberOfExternals();

			if ( maxOpenShells<minOpenShells )
				goto ende;

/*			for ( INT k=0 ; k<molist.getNumberOfInternals() ; k++ )
				cout << j[k] << " ";
			cout << endl;
*/

//			cout << conf[i+1] << endl;

		RepDiag<MatrixType, VectorType>	repDiag(repMatsP5, conf[i+1], minOpenShells, maxOpenShells);

		EnergyMap	energyMap(PTReference->getNumberOfRoots(), 
			mrinp.getEstimationMode());
		
//			,molist);

//			goto ende;

INT	selectAll = (molist.getNumberOfExternals()==0) && selectInternal;


			//	calculate excitation levels relative to references
			{
			const ConfigurationSet	*RefConfSet = &mrinp.getRefConfSet();
			INT	level = 10000000;
			
			Pix	iter = RefConfSet->first();
				while ( iter )
				{
				INT	l = Configuration<MOType>::calcExcitationOrder((*RefConfSet)(iter), conf[i+1]);
					
					if ( l<level )
						level = l;
					RefConfSet->next(iter);
				}
				selectAll |= !!(selectNthExcitation & (1<<(level-1)));
			}
			
//			cout << "!!!!!!!!!!!!!!!!!!!!!!" << endl;

			for ( 
				MRTreeIterator iter = PTReference->getPT0Wave()->firstInTree() ;
				!PTReference->getPT0Wave()->isLastInTree(iter) ;
				PTReference->getPT0Wave()->nextInTree(iter) 
				)
			{
/*				if ( molist.getNumberOfExternals()==0 )
				{
					dc.calcDiffConf(
						PTReference->getPT0Wave()->getConfigurationSAFNr(iter),
						conf[i+1]);
					tablecase->calcLess3(dc, 
							CoulombIndex, ExchangeIndex);
					HMatElements<MatrixType, VectorType>	*repMats = 
						(HMatElements<MatrixType, VectorType> *)
						(*cache)[TableKey(*tablecase)];

					RepMatType	dE[PTReference->getNumberOfRoots()*
						repMats->getNumberOfColumns()];

						repMats->PTcalc(dc, *tablecase, 
							CoulombIndex, ExchangeIndex, dE,
							*PTReference, iter.i, &repDiag);

				}
				else
*/				{
//				cout << "match" << endl;
				MatchingMOListIterator	match(
						PTReference->getPT0Wave()->getConfigurationSAFNr(iter),
						conf[i+1], &molist, 
						minOpenShells - conf[i+1].getNumberOfOpenShells());
	
//				cout << "match ok" << endl;
					
					while ( !match.isEnd() )
					{
//					cout << "jjjjjjjjjjjjjjjjjjjjj" << endl;
						dcGen = match.getGenDiffConf();
//					cout << dcGen << endl;
/*					DiffConf<GeneralizedMO>	dcGen1(dcGen);
						dcGen1.simplify();
					cout << "++++++++++++++++++++" << endl;
*/						tablecaseGen->calcLess3(dcGen, 
							CoulombIndexGen, ExchangeIndexGen);
						CoulombIndexT = CoulombIndexGen;
						ExchangeIndexT = ExchangeIndexGen;
//							cout << dcGen << endl;
//							cout << *tablecaseGen << endl << endl;
//							cout << CoulombIndexT << endl;
//							cout << ExchangeIndexT << endl;
						HMatElements<MatrixType, VectorType>	*repMats = 
						(HMatElements<MatrixType, VectorType> *)
							(*cache)[TableKey(*tablecaseGen)];

//						cout << "++++++++++++++++++++" << endl;
					MOListIterator	moiter = match.getMOListIterator();
//						cout << "++++++++++++++++++++" << endl;

//					RepMatType	dE[PTReference->getNumberOfRoots()*
//						repMats->getNumberOfColumns()*moiter.getN()];


//						cout << moiter << endl;
//						cout << match.getExtNotRunning() << endl;

//						cout << "AAAAAAAAA" << endl;
				
						repDiag.addExt(match.getExtNotRunning(), 
							dcGen.getSame().getNumberOfOpenShells() +
							dcGen.getTo().getNumberOfOpenShells(),
							dcGen.getPosTo(dcGen.getTo().getNthExtPosOpen(0)),
							dcGen.getPosTo(dcGen.getTo().getNthExtPosOpen(1))
							);

//						cout << "BBBBBBBBB" << endl;

/*
				Configuration<GeneralizedMO> cg = Configuration<GeneralizedMO>(dcGen.getSame());
						cg.create(dcGen.getTo());
				MOListIterator	it = moiter;
				
					if ( it.getN() )
					do {
//						cout << cg << endl;
					Configuration<MOType> tconf = Configuration<MOType>(cg, it);
//						cout << it << endl;
//						cout << tconf << endl;
					double	de = calcSingle(PTReference->getPT0Wave()->getConfigurationSAFNr(iter),
							tconf);
						cout << "~a " << setiosflags(ios::scientific) << de << endl;
					}
					while ( it.next() );
*/					
						
				
						repMats->PTcalc(dcGen, *tablecaseGen, 
							moiter, CoulombIndexT, ExchangeIndexT, 
							energyMap, *diagHist,
							match.getExtNotRunning(),
							*PTReference, iter.i, &repDiag);
//						cout << "CCCCCCCCC" << endl;
							
/*						for ( INT kk=0 ; kk <PTReference->getNumberOfRoots()*
						repMats->getNumberOfColumns()*moiter.getN() ; kk++ )
							cout << dE[kk] << " ";
						cout << endl;
*/
						match.next();
					}
//					cout << "OK1" << endl;
				}
//				cout << "OK2" << endl;
			}		

//*::::::::::::::::::::
			if ( preSelConfs )
			{
				energyMap.select(tupelWInt,
					preSelConfs, molist.getNumberOfExternals(),
					preESum);
			}
			else
				energyMap.select(tupelWInt,
					mrinp.getMORestriction(),
					mrinp.getMOStatistics(),
					enlargeRef,
					Threshold,
					mrinp.getRandThresh(),
					PTReference->getHistogramP(),
					PTReference->getPTTotalSum(),
					PTReference->getPTSumP(),
					genConfCount[molist.getNumberOfExternals()],
					genCSFCount[molist.getNumberOfExternals()],
					genRestrictedConfCount[molist.getNumberOfExternals()],
					genRestrictedCSFCount[molist.getNumberOfExternals()],
					RestrictedESum,
					
					selectAll);
//::::::::::::::::::::::*/
//			cout <<	"gen:" << energyMap.getNumberOfConfigurations();
//			cout << "sel:" << energyMap.select(1.0e-3) << endl;
		}
ende:
//		cout << "EEEEEEEEEEEEE" << endl;

		if ( molist.getNumberOfInternals()==0 )
			break;
	}
delete[] conf;
//	cout << "ENDE" << endl;
}


template <class MatrixType, class VectorType>
EnergyType	Selector<MatrixType, VectorType>::calcSingle(
	Configuration<MOType> mainA,
	Configuration<MOType> mainB)
{
	if ( Configuration<MOType>::calcExcitationOrder(mainA, mainB)>2 )
		return 0;

	dc.calcDiffConf(mainA, mainB);
	this->tablecase->calcLess3(dc);

	HMatElements<MatrixType, VectorType>	*repMats = (HMatElements<MatrixType, VectorType> *)
		(*cache)[TableKey(*this->tablecase)];


//	cout << "mainA=" << mainA << endl;
	cout <<  mainB << endl;
/*	cout << "diffConf=" << dc << endl;
	cout << "tablecase=" << *tablecase << endl;
	cout << "repMats:" << *repMats << endl;
*/
INT	r = repMats->getNumberOfRows();
INT	c = repMats->getNumberOfColumns();
EnergyType	*p = new EnergyType[r*c];

	repMats->getMatrix(p, dc, *this->tablecase);


	dc.calcDiffConf(mainB, mainB);
	this->tablecase->calcLess3(dc);
	repMats = (HMatElements<MatrixType, VectorType> *)
		(*cache)[TableKey(*this->tablecase)];
EnergyType	*pDiag = new EnergyType[c*c];

	repMats->getMatrix(pDiag, dc, *this->tablecase);


INT	root = 0;
INT	PTRefNr = 0;

EnergyType	e = 0;
	for ( INT i=0 ; i<c ; i++ )
	{
	EnergyType	ee = 0;
		for ( INT j=0 ; j<r ; j++ )
		{
//			cout << PTReference->getCoef(PTRefNr, root, j) << " ";
			ee += PTReference->getCoef(PTRefNr, root, j) * 
				p[j*c + i];
		}
//		cout << endl;
//		cout << "ee=" << ee << endl;
//		cout << "div=" << (PTReference->getRoot(root) - pDiag[(1+c)*i]) << endl;
		e += ee*ee / (PTReference->getRoot(root) - pDiag[(1+c)*i]);
	}

	delete p;
	delete pDiag;
	return e;
}



template <class MatrixType, class VectorType>
void	Selector<MatrixType, VectorType>::printResults() const
{
	cout << "=================================================================" << endl;
	cout << endl;
	cout << "                      R e s u l t s" << endl;	
	cout << endl;
	cout << "=================================================================" << endl;
	cout << endl;
	cout << endl;
	
	if ( preSelConfs )
	{
		for ( INT i=0 ; i<PTReference->getNumberOfRoots() ; i++ )
		{
			cout << "root no. " << i+1 << ": " << PTReference->getRoot(i) <<
				": " << preESum[i]*1e3 << " mH" << endl;
		}
	}
	else
	{
		cout << *PTReference->getPTTotalSum() << endl << endl;
		if ( mrinp.getMORestriction() )
		{
			cout << endl;
			cout << endl;
			cout << "energy sum from restricted configurations:" << endl;
			cout << endl;
			for ( INT i=0 ; i<PTReference->getNumberOfRoots() ; i++ )
				cout << "  root no. " << i+1 << "/mH";
			cout << endl;
			for ( INT i=0 ; i<PTReference->getNumberOfRoots() ; i++ )
				cout << "---------------";
			cout << endl;
			for ( INT i=0 ; i<PTReference->getNumberOfRoots() ; i++ )
				cout << setw(15) << RestrictedESum[i]*1000;
			cout << endl;
			cout << endl;
		}


		if ( verbosity.isActive(Verbosity::SelectionPerRoot) )
			for ( INT i=0 ; i<PTReference->getNumberOfRoots() ; i++ )
			{
				cout << "root no. " << i+1 << ": " << PTReference->getRoot(i) << endl;
				PTReference->getHistogram(i)->printPrettyASCII();
				cout << endl;
				cout << *PTReference->getPTSum(i) << endl;	
				cout << endl;
				cout << endl;
			}
	}

	if ( verbosity.isActive(Verbosity::DiagHist) )
	{
		cout << endl << endl << endl;
		cout << "diagonal histogram:" << endl;
		diagHist->printPrettyASCII();
		cout << endl << endl << endl;
	}
}


/*
template class Selector<float, float>;
template class Selector<double, float>;
template class Selector<float, double>;
*/
template class Selector<double, double>;
