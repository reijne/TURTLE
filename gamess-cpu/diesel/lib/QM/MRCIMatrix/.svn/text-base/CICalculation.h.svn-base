//***********************************************************************
//
//	Name:			CICalculation.h
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


#ifndef __CICalculation_h
#define __CICalculation_h

#include "../../../config.h"

#include <math.h>



#include "../RepresentationMatrices/HMatElements.h"

#include "../IntegralContainer/FourIndexIntegralContainer.h"
#include "../IntegralContainer/TwoIndexIntegralContainer.h"

#include "../IntegralContainer/RIFourIndexIntegralContainer.h"
#include "../IntegralContainer/RITwoIndexIntegralContainer.h"

#include "../IO/Fortran/Fort31File.h"

template <class MatrixType> class RepresentationMatrixFortranInterface;
template <class TMOType> class TableCase;

class SharedMemory;
class MRMOs;

template <class MatrixType, class VectorType>
class CICalculation {
public:
enum IntegralStorageMode { storeNone, storeAll };

	CICalculation(INT Multiplicity,
		INT numberOfSlaves, MRMOs mrmos, Fort31File fort31File,
		IntegralStorageMode);
	~CICalculation();


	INT	getMultiplicity() const 
	{	return repMatFInt->getMultiplicity();	}

	void	freeIntegrals();


protected:

BinomialCoefficient	*binom;					// binomial coefficients
TableCase<MOType>	*tablecase;				// table case



FourIndexIntegralContainer	*intContainer;	// 4 index integrals
TwoIndexIntegralContainer	*int2Container;	// 2 index integrals

SharedMemory	*Int4sharedMem;				// shared memory containing 
											// 4 index integrals


Fort31File	Fort31;							// file from which the integrals are
											// read

RepresentationMatrixFortranInterface<MatrixType>	*repMatFInt;	
										// interface to FORTRAN calculation
										// of representation matrices


MRMOs mrmos;

INT	freeOnDestruct;						// flag if destructors on contained
										// objects should be called 
};



#endif
