//***********************************************************************
//
//	Name:			EnlargeReferenceSpace.cc
//
//	Description:	
//
//	Author:			Michael Hanrath
//
//	Version:		0.0
//
//	Date:			01.04.1998
//
//
//
//
//
//***********************************************************************

#include "EnlargeReferenceSpace.h"

#include "../../../lib/Math/FortranNumeric/FortranEigenProblems.h"
#include "../../../lib/QM/MRTree/Set/NExternalsSet.h"
#include "../../../lib/QM/Cache/VarSizeReadOnlyCache.h"
#include "../../../lib/QM/MRTree/MRTreeIterator.h"
#include "MRConfInput.h"
#include "../../../lib/QM/Math/SpinEigenFunctionDegeneration.h"

#include <iostream>
#include <iomanip>

using namespace std;

template <class RepMatType, class VectorType>
EnlargeReferenceSpace<RepMatType, VectorType>::EnlargeReferenceSpace(
	const MRConfInput &mrinp)
{
	PT0Wave = new NExternalsSet(
		mrinp.getMRMOs(),
		mrinp.getNumberOfElectrons(),
		mrinp.getMultiplicity(),
		mrinp.getNumberOfRoots(),
		((MRConfInput) mrinp).getPTRefConfSet()
		);


		
	dim = PT0Wave->getNumberOfTotalSpinAdaptedFunctions();


	pRefMat = new RepMatType[dim*(dim+1)/2];

DiffConf<MOType>	diffConf;


	spinEigs = new SpinEigenFunctionDegeneration(mrinp.getMultiplicity(), MAXOPENSHELLS);
	binom = new BinomialCoefficient(MAXOPENSHELLS);
	tablecase = new TableCase<MOType>(binom);
	cache = new VarSizeReadOnlyCache<TableKey, HMatElements<RepMatType, VectorType> >
		(2048, (1<<20));
	

	memset(pRefMat, 0, dim*(dim+1)/2*sizeof(RepMatType));
	
	for ( MRTreeIterator iteri = PT0Wave->firstInTree() ;
		!PT0Wave->isLastInTree(iteri) ; PT0Wave->nextInTree(iteri) )
	{
	ConfigurationSAFNr<MOType> mainA = PT0Wave->getConfigurationSAFNr(iteri);
	INT	safsA = mainA.getSAFInc();

		for ( MRTreeIterator iterj = PT0Wave->firstInTree() ;
			iterj.i<=iteri.i ; PT0Wave->nextInTree(iterj) )
		{
		ConfigurationSAFNr<MOType> mainB = PT0Wave->getConfigurationSAFNr(iterj);

			if ( ConfigurationSAFNr<MOType>::calcExcitationOrder(mainA, mainB)>2 )
				continue;
				
			diffConf.calcDiffConf(mainA, mainB);
			tablecase->calcLess3(diffConf);

			HMatElements<RepMatType, VectorType>	*repMats = (HMatElements<RepMatType, VectorType> *)
				(*cache)[TableKey(*tablecase)];


		INT	safsB = mainB.getSAFInc();
		
		RepMatType	*p = new RepMatType[safsA*safsB];
		
			repMats->getMatrix(p, diffConf, *tablecase);

		RepMatType	*pp = p;
			if ( iteri.i == iterj.i )
			{
				for ( INT k=0 ; k<safsA ; k++ )
				{
					for ( INT l=0 ; l<=k ; l++ )
						pRefMat[iover2(mainA.getSAFNr()+k)
						+ mainB.getSAFNr()+l] = *pp++;
					pp += safsA - (k + 1);
				}
			}
			else
			{
				for ( INT k=0 ; k<safsA ; k++ )
					for ( INT l=0 ; l<safsB ; l++ )
						pRefMat[iover2(mainA.getSAFNr()+k)
						+ mainB.getSAFNr()+l] = *pp++;
			}		
			delete p;
		}
	}
	


RepMatType	*lambda = new RepMatType[dim];
RepMatType	*alpha = new RepMatType[dim*dim];

RepMatType	*pRefMatCopy = new RepMatType[dim*(dim+1)/2];
	memcpy(pRefMatCopy, pRefMat, dim*(dim+1)/2*sizeof(RepMatType));


	hqrii1(
		&dim, &dim, pRefMatCopy, lambda, alpha, &dim
	);
	


	// store roots
	nRoots = mrinp.getNumberOfRoots();
	refEnergies = new RootType[nRoots];
	roots = new INT[nRoots];
	for ( INT i=0, j=0 ; i<dim && j<nRoots ; i++ )
		if  ( mrinp.getRootNumber(j)-1==i )
			roots[j++] = i;

	for ( INT i=0 ; i<nRoots ; i++ )
		refEnergies[i] = lambda[mrinp.getRootNumber(i)-1];


	delete lambda;
	delete alpha;


	enlargedEnergies = new EnergyType[nRoots];
	memset(enlargedEnergies, 0, nRoots*sizeof(EnergyType));
}


template <class RepMatType, class VectorType>
EnlargeReferenceSpace<RepMatType, VectorType>::~EnlargeReferenceSpace()
{
	delete refEnergies;
	delete roots;
	delete enlargedEnergies;
	delete pRefMat;
	delete PT0Wave;

	delete cache;
	delete tablecase;
	delete binom;
	delete spinEigs;
}



template <class RepMatType, class VectorType>
const EnergyType *	EnlargeReferenceSpace<RepMatType, VectorType>::calcEnlargementEnergy(
	const Configuration<MOType> &main)
{
INT	safs = (*spinEigs)(main.getNumberOfOpenShells());
INT	dim2 = dim + safs;


	memset(enlargedEnergies, 0, nRoots*sizeof(EnergyType));

RepMatType	*lambda = new RepMatType[dim2];
RepMatType	*alpha = new RepMatType[dim2*dim2];

RepMatType	*pRefMatCopy = new RepMatType[dim2*(dim2+1)/2];
	memset(pRefMatCopy, 0, dim2*(dim2+1)/2*sizeof(RepMatType));
	memcpy(pRefMatCopy, pRefMat, dim*(dim+1)/2*sizeof(RepMatType));


DiffConf<MOType>	diffConf;

	for ( MRTreeIterator iteri = PT0Wave->firstInTree() ;
		!PT0Wave->isLastInTree(iteri) ; PT0Wave->nextInTree(iteri) )
	{
	ConfigurationSAFNr<MOType> mainA = PT0Wave->getConfigurationSAFNr(iteri);
	INT	safsA = mainA.getSAFInc();

		if ( ConfigurationSAFNr<MOType>::calcExcitationOrder(mainA, main)>2 )
			continue;

		diffConf.calcDiffConf(mainA, main);
		tablecase->calcLess3(diffConf);
		if ( tablecase->getP()==5 )
			return enlargedEnergies;

		HMatElements<RepMatType, VectorType>	*repMats = (HMatElements<RepMatType, VectorType> *)
			(*cache)[TableKey(*tablecase)];



	RepMatType	*p = new RepMatType[safsA*safs];

		repMats->getMatrix(p, diffConf, *tablecase);

	RepMatType	*pp = p;
		for ( INT k=0 ; k<safs ; k++ )
			for ( INT l=0 ; l<safsA ; l++ )
			{
//				cout << k << " " << l << " " << iover2(dim+k) << " " << mainA.getSAFNr() << " " << *pp << endl;
				pRefMatCopy[iover2(dim+k) + mainA.getSAFNr() + l] = *pp++;
			}
		delete p;
	}

	{
		diffConf.calcDiffConf(main, main);
		tablecase->calcLess3(diffConf);

		HMatElements<RepMatType, VectorType>	*repMats = (HMatElements<RepMatType, VectorType> *)
			(*cache)[TableKey(*tablecase)];

	RepMatType	*p = new RepMatType[safs*safs];

		repMats->getMatrix(p, diffConf, *tablecase);

	RepMatType	*pp = p;
		for ( INT k=0 ; k<safs ; k++ )
		{
			for ( INT l=0 ; l<=k ; l++ )
				pRefMatCopy[iover2(dim+k) + dim + l] = *pp++;
			pp += safs - (k + 1);
		}
	}
	
/*	cout << "!!!!!!!!!!!!!!!!!!reference matrix!!!!!!!!!!!!!!!!!!!!:" << endl;
	for ( INT i=0 ; i<dim2 ; i++ )
	{
		for ( INT j=0 ; j<dim2 ; j++ )
			if ( j<=i )
				cout << setw(14) << setprecision(6) << pRefMatCopy[i*(i+1)/2 + j] << " ";
			else
				cout << "\t";
		cout << endl;
	}
*/

	
/*
RepMatType	*ppp = new RepMatType[dim2*(dim2+1)/2];
	memcpy(ppp, pRefMatCopy, dim2*(dim2+1)/2*sizeof(RepMatType));

*/
	hqrii1(
		&dim2, &dim2, pRefMatCopy, lambda, alpha, &dim2
	);
	
	for ( INT i=0 ; i<nRoots ; i++ )
		enlargedEnergies[i] = lambda[roots[i]] - refEnergies[i];


/*	if ( enlargedEnergies[0]==0 )
	{
		cout << main << endl;
		cout << "!!!!!!!!!!!!!!!!!!reference matrix!!!!!!!!!!!!!!!!!!!!:" << endl;
		for ( INT i=0 ; i<dim2 ; i++ )
		{
			for ( INT j=0 ; j<dim2 ; j++ )
				if ( j<=i )
					cout << setw(14) << setprecision(6) << ppp[i*(i+1)/2 + j] << " ";
				else
					cout << "\t";
			cout << endl;
		}
		cout << endl;
		cout << endl;
		cout << endl;
	}
*/


	delete lambda;
	delete alpha;
	delete pRefMatCopy;
	return enlargedEnergies;
}



template <class RepMatType, class VectorType>
ostream & operator << (ostream &s, const EnlargeReferenceSpace<RepMatType, VectorType> &e)
{
	s << "reference matrix:" << endl;
	for ( INT i=0 ; i<e.dim ; i++ )
	{
		for ( INT j=0 ; j<e.dim ; j++ )
			if ( j<=i )
				s << setw(14) << setprecision(6) << e.pRefMat[i*(i+1)/2 + j] << " ";
			else
				s << "\t";
		s << endl;
	}
	s << endl;
	s << "eigenvalues:" << endl;
	for ( INT i=0 ; i<e.nRoots ; i++ )
			s << setw(2) << e.roots[i] 
				<< setprecision(6) << setw(14) << e.refEnergies[i] << endl; 
	return s;
}



template class EnlargeReferenceSpace<float, float>;
template class EnlargeReferenceSpace<double, float>;
template class EnlargeReferenceSpace<float, double>;
template class EnlargeReferenceSpace<double, double>;
