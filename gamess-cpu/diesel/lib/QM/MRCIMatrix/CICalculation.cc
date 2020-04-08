//***********************************************************************
//
//	Name:			CICalculation.cc
//
//	Description:	
//
//	Author:			Michael Hanrath
//
//	Version:		0.0
//
//	Date:			09.09.1998
//
//
//
//
//***********************************************************************

#include "CICalculation.h"

#include "../Parallel/SharedMemory.h"
#include "../MO/MRMOs.h"
#include "../IO/Verbosity.h"

#include "../IO/Fortran/Fort31FirstRecord.h"


template <class MatrixType, class VectorType>
CICalculation<MatrixType, VectorType>::CICalculation(
	INT Multiplicity,
	INT numberOfSlaves, const MRMOs _mrmos, Fort31File _Fort31,
	IntegralStorageMode integralStorageMode)
{
	freeOnDestruct = 1;
	Fort31 = _Fort31;
	mrmos = _mrmos;
	
	binom = new BinomialCoefficient(MAXOPENSHELLS);
	tablecase = new TableCase<MOType>(binom);

	repMatFInt = 
		new RepresentationMatrixFortranInterface<MatrixType>(Multiplicity);




	switch ( integralStorageMode ) {
	case storeAll:
		{
//----------------------------------------------------------------------
//	load integrals
	
//----------------------
//  take correct Container for RI-Calculation
//---------------------- 
			if(Fort31.format != RIFormat)
				intContainer = new FourIndexIntegralContainer(mrmos);
			else
			{
				intContainer = new RIFourIndexIntegralContainer(mrmos,100);
				//intContainer = new RIFourIndexIntegralContainer(mrmos);
				//intContainer->precalcMode=storeAll;
			}
		INT	size = intContainer->getDataSize();

			if ( numberOfSlaves==0 )
				Int4sharedMem = new SharedMemory(size, SharedMemory::StandAlone);
			else
				Int4sharedMem = new SharedMemory(size, SharedMemory::Master);

			intContainer->setSharedMem(Int4sharedMem);


			if ( verbosity.isActive(Verbosity::Integrals) )
			{
				cout << "total number of two electron integral containers: " << 
					intContainer->getNumberOfLeaves() << endl;
				cout << "total number of two electron integrals: " << 
					intContainer->getNumberOfIntegrals() << endl;

				cout << "allocating memory for two electron integrals..." << std::flush;
			}

			if ( intContainer->allocate() )
			{
				if ( verbosity.isActive(Verbosity::Integrals) )
					cout << "success" << endl;
			}
			else
			{	cout << "failed" << endl;
				exit(1);
			}

			intContainer->loadIntegrals(Fort31);

//--------------------------------------------------------------------------

//----------------------
//  take correct Container for RI-Calculation
//----------------------
		if(Fort31.format != RIFormat)
			int2Container = new TwoIndexIntegralContainer(mrmos);
		else
			int2Container = new RITwoIndexIntegralContainer(mrmos);

		double	core = int2Container->loadIntegrals(Fort31);


//--------------------------------------------------------------------------

		Fort31FirstRecord	f31(Fort31);


			cout << "VNUC = " << f31.getVNUC() << endl;
			cout << "ZERO = " << f31.getZERO() << endl;
			cout << "CORE = " << core << endl;



//--------------------------------------------------------------------------
//	initialize static data in HMatElements<MatrixType, VectorType>
		HMatElements<MatrixType, VectorType> init1(intContainer, int2Container, core);
		RepresentationMatrix<MatrixType> init2;
		}
		break;
		
	case storeNone:
		{
			intContainer = NULL;
			int2Container = NULL;
			Int4sharedMem = NULL;

//--------------------------------------------------------------------------
//	initialize static data in HMatElements<MatrixType, VectorType>
		HMatElements<MatrixType, VectorType> init1(mrmos.getMaxMO());
		RepresentationMatrix<MatrixType> init2;
		}
		break;
	}
	if ( verbosity.isActive(Verbosity::MOs) )
		cout << endl << mrmos << endl;
}



template <class MatrixType, class VectorType>
CICalculation<MatrixType, VectorType>::~CICalculation()
{
	if ( freeOnDestruct )
	{
		if ( intContainer )
			delete intContainer;

		if ( int2Container )
			delete int2Container;

		delete repMatFInt;
		delete tablecase;
		delete binom;


		if ( Int4sharedMem )
			delete Int4sharedMem;
	}
}


template <class MatrixType, class VectorType>
void	CICalculation<MatrixType, VectorType>::freeIntegrals()
{
	if ( Int4sharedMem )
	{
		delete Int4sharedMem;
		Int4sharedMem = NULL;
	}

	if ( intContainer )
	{
		delete intContainer;
		intContainer = NULL;
	}
}


template class CICalculation<float, float>;
template class CICalculation<double, float>;
template class CICalculation<float, double>;
template class CICalculation<double, double>;
