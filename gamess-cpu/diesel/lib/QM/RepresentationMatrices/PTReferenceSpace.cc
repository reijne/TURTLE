//***********************************************************************
//
//	Name:			PTReferenceSpace.cc
//
//	Description:	
//
//	Author:			Michael Hanrath
//
//	Version:		0.0
//
//	Date:			08.04.1997
//
//
//
//
//
//***********************************************************************

#include "PTReferenceSpace.h"


#include "RepresentationMatrices.h"
#include "../IO/TimeTicks.h"
#include "../MRTree/MRTreeIterator.h"
#include "../Configuration/Configuration.h"
#include "../Cache/VarSizeReadOnlyCache.h"
#include "../../Math/etc/iover2.h"
#include "../MO/Iterators/MOIterator.h"

#include "MainRefSet.h"

#include "../../Math/FortranNumeric/FortranEigenProblems.h"

#include "../MRTree/Set/NExternalsSet.h"
#include "../MRTree/EnergyMap/PTSum.h"
#include "../MRCIMatrix/MRCIMatrix.h"
#include "../MRTree/Diag/NExternalsDiag.h"
#include "../IntegralContainer/MRFockMatrix.h"
#include "../MRCIMatrix/MRMPH0Matrix.h"
#include "../Davidson/Roots.h"
#include "../Davidson/ConjGrad.h"
#include "../Davidson/SLEJacobi.h"
#include "../../Container/DiskBuffer.h"


#include "../Davidson/CIMat/CIEig.h"
#include "../Davidson/CIMat/CIMatFull.h"
#include "../Davidson/CIMat/CIMatSparse.h"

#include "../IO/Verbosity.h"

#include <iomanip>
#include <fstream>

template <class MatrixType, class VectorType>
PTReferenceSpace<MatrixType, VectorType>::PTReferenceSpace(MRConfInput mrinp,
	const CICalculation<MatrixType, VectorType> & ciCalculation)
{
	PT0Wave = new NExternalsSet(
		mrinp.getMRMOs(),
		mrinp.getNumberOfElectrons(),
		mrinp.getMultiplicity(),
		mrinp.getNumberOfRoots(),
		mrinp.getPTRefConfSet()
		);
	cout << endl;

/*	for ( MRTreeIterator i = PT0Wave->firstInTree() ;
			!PT0Wave->isLastInTree(i) ; PT0Wave->nextInTree(i))
	{
		cout << i.Citer[0].i<< " " << i.Citer[1].i<< " " << i.Citer[2].i<< " " << i.Citer[3].i << ": ";
		cout << PT0Wave->getConfigurationSAFNr(i) << endl;
		cout << PT0Wave->getConfigurationSAFNr(i).getNExt() << endl << endl;
	}
*/
	cout << "zero order wave function ci-problem (basis of perturbation theory):"
		<< endl;
	cout << "number of configurations: " << PT0Wave->getNumberOfLeaves() << endl;
	cout << "dimension of ci matrix  : " 
		<< PT0Wave->getNumberOfTotalSpinAdaptedFunctions() << endl;
		
INT	dim = PT0Wave->getNumberOfTotalSpinAdaptedFunctions();


	if ( mrinp.getNumberOfRoots()>dim )
	{
		cout << "error: too many (" << mrinp.getNumberOfRoots()
			<< ") roots ordered" << endl;
		exit(1);
	}


MatrixType	*pRefMat = new MatrixType[dim*(dim+1)/2];

DiffConf<MOType>	diffConf;
TimeTicks	ticks;




BinomialCoefficient	*binom = new BinomialCoefficient(MAXOPENSHELLS);
TableCase<MOType>	*tablecase = new TableCase<MOType>(binom);
VarSizeReadOnlyCache<TableKey, HMatElements<double, double> >	*cache = 
		new VarSizeReadOnlyCache<TableKey, HMatElements<double, double> >
		(2048, (1<<20));
	
	cout << endl << "generating ci-matrix... " << flush;
	//	not very efficient
	ticks.start();
	memset(pRefMat, 0, dim*(dim+1)/2*sizeof(MatrixType));
	
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

			HMatElements<MatrixType, VectorType>	*repMats = 
				(HMatElements<MatrixType, VectorType> *)
				(*cache)[TableKey(*tablecase)];


		INT	safsB = mainB.getSAFInc();
		
		MatrixType	*p = new MatrixType[safsA*safsB];
		
			repMats->getMatrix(p, diffConf, *tablecase);

		MatrixType	*pp = p;
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
			delete[] p;
		}
	}
	ticks.stop();
	cout << ticks << endl;
	
	delete cache;
	delete tablecase;
	delete binom;

	
	if ( verbosity.isActive(Verbosity::RefMat) )
	{
		cout << "reference matrix:" << endl;
		for ( INT i=0 ; i<dim ; i++ )
		{
			for ( INT j=0 ; j<dim ; j++ )
				if ( j<=i )
					cout << setw(14) << pRefMat[i*(i+1)/2 + j] << " ";
				else
					cout << "\t";
			cout << endl;
		}
		cout << endl;
	}


	cout.precision(8);
	cout.setf(ios::fixed);
MatrixType	*lambda = new MatrixType[dim];
MatrixType	*alpha = new MatrixType[dim*dim];

	cout << endl << "diagonalizing ci-matrix... " << flush;
	ticks.start();

	hqrii1(
		&dim, &dim, pRefMat, lambda, alpha, &dim
	);
	
	ticks.stop();
	cout << ticks << endl;


	if ( verbosity.isActive(Verbosity::RefMatEigenValues) )
	{
		cout << "eigenvalues of reference matrix:" << endl;
		for ( INT i=0, j=0 ; i<dim ; i++ )
		{
			cout << lambda[i];
			if  ( j<mrinp.getNumberOfRoots() && mrinp.getRootNumber(j)-1==i )
			{
				j++;
				cout << " <--- selected";
			}
			cout << endl;

		}
		cout << endl;
	}

	if ( verbosity.isActive(Verbosity::RefMatEigenVectors) )
	{
	INT	 nPrint = (INT) ceil(mrinp.getNumberOfRoots()*1.5);
		if ( nPrint>dim )
			nPrint = dim;
		cout << "eigenvectors of reference matrix (columnwise):" << endl;
		cout << setw(5) << "ref #  ";
		for ( INT j=0 ; j<nPrint ; j++ )
			cout << setw(6) << "root #" << setw(3) << j+1 << "  ";
		cout << endl;
	INT	safc = 1;
	MRTreeIterator iteri = PT0Wave->firstInTree();
	INT	k = 0;
		for ( INT i=0 ; i<dim ; i++ )
		{
			if ( !--safc )
			{
				safc = PT0Wave->getConfigurationSAFNr(iteri).getSAFInc();
				PT0Wave->nextInTree(iteri);
				cout << setw(5) << ++k << " ";
			}
			else
				cout << setw(5) << "|" << " ";

			for ( INT j=0 ; j<nPrint ; j++ )
				cout << setprecision(5) << setw(10) << alpha[j*dim + i] << " ";
			cout << endl;
		}
		cout << endl;
	}
	
	
	// store roots
	nRoots = mrinp.getNumberOfRoots();
	roots = new RootType[nRoots];
	for ( INT i=0 ; i<nRoots ; i++ )
		if ( mrinp.getRootEnergyP() )
			roots[i] = mrinp.getRootEnergy(i);
		else
			roots[i] = lambda[mrinp.getRootNumber(i)-1];
	delete[] lambda;



	// store eigenvectors
MainRefSet	mains(nRoots, mrinp.getActiveReferenceThreshold());
//	refSet = mains.getMainRefSet();
	refSet = ConfigurationSet();
	coefs = new CoefType**[PT0Wave->getNumberOfLeaves()];
double* tsum = new double[nRoots];
	memset(tsum, 0, nRoots*sizeof(double));
	for ( MRTreeIterator iter = PT0Wave->firstInTree() ;
		!PT0Wave->isLastInTree(iter) ; PT0Wave->nextInTree(iter) )
	{
	ConfigurationSAFNr<MOType> main = PT0Wave->getConfigurationSAFNr(iter);
	INT	safs = main.getSAFInc();
		
		coefs[iter.i] = new CoefType*[nRoots];
	double* sum = new double[nRoots];
		memset(sum, 0, nRoots*sizeof(double));
		for ( INT i=0 ; i<nRoots ; i++ )
		{
			coefs[iter.i][i] = new CoefType[safs];
			for ( INT j=0 ; j<safs ; j++ )
			{
				coefs[iter.i][i][j] = 
					alpha[(mrinp.getRootNumber(i)-1)*dim + main.getSAFNr() + j];
				sum[i] += coefs[iter.i][i][j]*coefs[iter.i][i][j];

//				cout << mrinp.getRootNumber(i) << " " << main.getSAFNr() << " " << j << endl;
//				cout << "??????" << iter.i << " " << i << " " << j << " " << coefs[iter.i][i][j] << endl;
			}
			if ( sum[i]>mrinp.getActiveReferenceThreshold() )
			{
				refSet.add(main);
				tsum[i] += sum[i];
			}
		}
		mains.add(main, sum);
		delete[] sum;
	}

	cout << endl;	
	cout << endl;	
	for ( INT i=0 ; i<nRoots ; i++ )
		cout << "root " << i+1 << ": ci^2=" << tsum[i] << endl;
	cout << endl;	
	cout << endl;	
	cout << endl;	
	
//	mains.calc();
//	refSet = mains.getMainRefSet();
		

	mrFockMatrix = NULL;

	delete[] alpha;
	delete[] tsum;
	
	hist = new Histogram<EnergyType> * [nRoots];
	ptSum = new PTSum * [nRoots];
	for ( INT i=0 ; i<nRoots ; i++ )
	{
		hist[i] = new Histogram<EnergyType>(
			1e-16, 1e-2, Histogram<EnergyType>::Logarithmic, 40);
		ptSum[i] = new PTSum(
			mrinp.getNumberOfSelectionThresholds(), mrinp.getSelectionThresholdP(),
			mrinp.getEstimationMode()==EnergyMap::EpsteinNesbet
			);
	}
	ptTotalSum = new PTSum(
		mrinp.getNumberOfSelectionThresholds(), mrinp.getSelectionThresholdP(),
		mrinp.getEstimationMode()==EnergyMap::EpsteinNesbet,
		nRoots);
}

/*template <class MatrixType, class VectorType>
PTReferenceSpace<MatrixType, VectorType>::PTReferenceSpace(MRConfInput mrinp,
	const CICalculation<MatrixType, VectorType> & ciCalculation)
{
	PT0Wave = new NExternalsSet(
		mrinp.getMRMOs(),
		mrinp.getNumberOfElectrons(),
		mrinp.getMultiplicity(),
		mrinp.getNumberOfRoots(),
		mrinp.getPTRefConfSet()
		);
	cout << endl;


	cout << "zero order wave function ci-problem (basis of perturbation theory):"
		<< endl;
	cout << "number of configurations: " << PT0Wave->getNumberOfLeaves() << endl;
	cout << "dimension of ci matrix  : " 
		<< PT0Wave->getNumberOfTotalSpinAdaptedFunctions() << endl;
		
INT	dim = PT0Wave->getNumberOfTotalSpinAdaptedFunctions();


	if ( mrinp.getNumberOfRoots()>dim )
	{
		cout << "error: too many (" << mrinp.getNumberOfRoots()
			<< ") roots ordered" << endl;
		exit(1);
	}



DiffConf<MOType>	diffConf;
TimeTicks	ticks;




BinomialCoefficient	*binom = new BinomialCoefficient(MAXOPENSHELLS);
TableCase<MOType>	*tablecase = new TableCase<MOType>(binom);
VarSizeReadOnlyCache<TableKey, HMatElements<double, double> >	*cache = 
		new VarSizeReadOnlyCache<TableKey, HMatElements<double, double> >
		(2048, (1<<20));
	
	cout << endl << "generating ci-matrix... " << flush;
	//	not very efficient
	ticks.start();

INT	 nRef = (INT) ceil(mrinp.getNumberOfRoots()*2);
	if ( nRef>dim )
		nRef = dim;

INT	 nCalcRoots = (INT) ceil(mrinp.getNumberOfRoots()*1.5);
	if ( nCalcRoots>dim )
		nCalcRoots = dim;

//CIMatFull	cimat(nRef, dim);
CIMatSparse	cimat(nRef, dim);

	
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

			HMatElements<MatrixType, VectorType>	*repMats = 
				(HMatElements<MatrixType, VectorType> *)
				(*cache)[TableKey(*tablecase)];


		INT	safsB = mainB.getSAFInc();
		
		MatrixType	*p = new MatrixType[safsA*safsB];
		
			repMats->getMatrix(p, diffConf, *tablecase);

		MatrixType	*pp = p;
		INT	ASAF = mainA.getSAFNr();
		INT	BSAF = mainB.getSAFNr();
			if ( iteri == iterj )
			{
				for ( INT k=0 ; k<safsA ; k++ )
				{
					for ( INT l=0 ; l<=k ; l++ )
					{
						if ( fabs(*pp)>1e-20 )
							cimat.append(ASAF+k, BSAF+l, *pp);
						pp++;
					}
					pp += safsA - (k + 1);
				}
			}
			else
			{
				for ( INT k=0 ; k<safsA ; k++ )
					for ( INT l=0 ; l<safsB ; l++ )
					{
						if ( fabs(*pp)>1e-20 )
							cimat.append(ASAF+k, BSAF+l, *pp);
						pp++;
					}
			}		
			delete p;
		}
	}
	ticks.stop();
	cout << ticks << endl;
	
	delete cache;
	delete tablecase;
	delete binom;

	cimat.chooseRefMat();
	
	if ( verbosity.isActive(Verbosity::RefMat) )
	{
		cout << "reference matrix:" << endl;
		for ( INT i=0 ; i<dim ; i++ )
		{
			for ( INT j=0 ; j<dim ; j++ )
				if ( j<=i )
					cout << setw(14) << cimat(i, j) << " ";
				else
					cout << "\t";
			cout << endl;
		}
		cout << endl;
	}


	cout.precision(8);
	cout.setf(ios::fixed);

	cout << endl << "diagonalizing ci-matrix... " << flush;
	ticks.start();


CIEig	eig(cimat, nCalcRoots);

//	hqrii1(
//		&dim, &dim, pRefMat, lambda, alpha, &dim
//	);
	
	ticks.stop();
	cout << ticks << endl;


	if ( verbosity.isActive(Verbosity::RefMatEigenValues) )
	{
		cout << "eigenvalues of reference matrix:" << endl;
		for ( INT i=0, j=0 ; i<nCalcRoots ; i++ )
		{
			cout << eig[i];
			if  ( mrinp.getRootNumber(j)-1==i )
			{
				j++;
				cout << " <--- selected";
			}
			cout << endl;

		}
		cout << endl;
	}


	if ( verbosity.isActive(Verbosity::RefMatEigenVectors) )
	{
		cout << "eigenvectors of reference matrix (columnwise):" << endl;
		cout << setw(5) << "ref #  ";
		for ( INT j=0 ; j<nCalcRoots ; j++ )
			cout << setw(6) << "root #" << setw(3) << j+1 << "  ";
		cout << endl;
	INT	safc = 1;
	MRTreeIterator iteri = PT0Wave->firstInTree();
	INT	k = 0;
		for ( INT i=0 ; i<dim ; i++ )
		{
			if ( !--safc )
			{
				safc = PT0Wave->getConfigurationSAFNr(iteri).getSAFInc();
				PT0Wave->nextInTree(iteri);
				cout << setw(5) << ++k << " ";
			}
			else
				cout << setw(5) << "|" << " ";

			for ( INT j=0 ; j<nCalcRoots ; j++ )
				cout << setprecision(5) << setw(10) << eig.getEV(j, i) << " ";
			cout << endl;
		}
		cout << endl;
	}
	
	
	// store roots
	nRoots = mrinp.getNumberOfRoots();
	roots = new RootType[nRoots];
	for ( INT i=0 ; i<nRoots ; i++ )
		if ( mrinp.getRootEnergyP() )
			roots[i] = mrinp.getRootEnergy(i);
		else
			roots[i] = eig[mrinp.getRootNumber(i)-1];



	// store eigenvectors
MainRefSet	mains(nRoots, mrinp.getActiveReferenceThreshold());
//	refSet = mains.getMainRefSet();
	refSet = ConfigurationSet();
	coefs = new CoefType**[PT0Wave->getNumberOfLeaves()];
double	tsum[nRoots];
	memset(tsum, 0, nRoots*sizeof(double));
	for ( MRTreeIterator iter = PT0Wave->firstInTree() ;
		!PT0Wave->isLastInTree(iter) ; PT0Wave->nextInTree(iter) )
	{
	ConfigurationSAFNr<MOType> main = PT0Wave->getConfigurationSAFNr(iter);
	INT	safs = main.getSAFInc();
		
		coefs[iter.i] = new CoefType*[nRoots];
	double	sum[nRoots];
		memset(sum, 0, nRoots*sizeof(double));
		for ( INT i=0 ; i<nRoots ; i++ )
		{
			coefs[iter.i][i] = new CoefType[safs];
			for ( INT j=0 ; j<safs ; j++ )
			{
				coefs[iter.i][i][j] = 
					eig.getEV(mrinp.getRootNumber(i)-1, main.getSAFNr() + j);
				sum[i] += coefs[iter.i][i][j]*coefs[iter.i][i][j];

//				cout << mrinp.getRootNumber(i) << " " << main.getSAFNr() << " " << j << endl;
//				cout << "??????" << iter.i << " " << i << " " << j << " " << coefs[iter.i][i][j] << endl;
			}
			if ( sum[i]>mrinp.getActiveReferenceThreshold() )
			{
				refSet.add(main);
				tsum[i] += sum[i];
			}
		}
		mains.add(main, sum);
	}

	cout << endl;	
	cout << endl;	
	for ( INT i=0 ; i<nRoots ; i++ )
		cout << "root " << i+1 << ": ci^2=" << tsum[i] << endl;
	cout << endl;	
	cout << endl;	
	cout << endl;	
	
//	mains.calc();
//	refSet = mains.getMainRefSet();
		

	mrFockMatrix = NULL;

	
	
	hist = new Histogram<EnergyType> * [nRoots];
	ptSum = new PTSum * [nRoots];
	for ( INT i=0 ; i<nRoots ; i++ )
	{
		hist[i] = new Histogram<EnergyType>(
			1e-16, 1e-2, Histogram<EnergyType>::Logarithmic, 40);
		ptSum[i] = new PTSum(
			mrinp.getNumberOfSelectionThresholds(), mrinp.getSelectionThresholdP(),
			mrinp.getEstimationMode()==EnergyMap::EpsteinNesbet
			);
	}
	ptTotalSum = new PTSum(
		mrinp.getNumberOfSelectionThresholds(), mrinp.getSelectionThresholdP(),
		mrinp.getEstimationMode()==EnergyMap::EpsteinNesbet,
		nRoots);
}
*/



template <class MatrixType, class VectorType>
PTReferenceSpace<MatrixType, VectorType>::~PTReferenceSpace()
{
	if( mrFockMatrix )
		delete mrFockMatrix;
	for ( INT i=0 ; i<PT0Wave->getNumberOfLeaves() ; i++ )
	{
		for ( INT j=0 ; j<nRoots ; j++ )
			delete coefs[i][j];
		delete coefs[i];
	}
	delete coefs;
	delete roots;
	delete PT0Wave;
	delete ptTotalSum;
	for ( INT i=0 ; i<nRoots ; i++ )
	{
		delete hist[i];
		delete ptSum[i];
	}
	delete hist;
	delete ptSum;
}



//template class PTReferenceSpace<float, float>;
//template class PTReferenceSpace<double, float>;
//template class PTReferenceSpace<float, double>;
template class PTReferenceSpace<double, double>;



