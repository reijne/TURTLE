//***********************************************************************
//
//	Name:			RepresentationMatrices.h
//
//	Description:	stores and handles the 
//					representation matrices of the 
//					symmetric group approach (SGA)
//
//					
//
//	Author:			Michael Hanrath
//
//	Version:		0.0
//
//	Date:			23.08.1996
//
//
//
//
//
//***********************************************************************

#ifndef __RepresentationMatrices_h
#define __RepresentationMatrices_h

#include "../../../config.h"

#include "CIVectors.h"

#include "../Cache/CacheableObject.h"
#include "../Cache/TableKey.h"

#include "../MO/MOType.h"
#include "../Configuration/TableCase.h"

#include "../IntegralContainer/IntegralType.h"

#include "../Configuration/Configuration.h"
#include "../Configuration/DiffConf.h"

#include "RepresentationMatrix.h"
#include "RepresentationMatrixFortranInterface.h"

#include "../IntegralContainer/FourIndexIntegralContainer.h"
#include "../IntegralContainer/TwoIndexIntegralContainer.h"

#include "../MO/Iterators/MOListIterator.h"
#include "../IntegralIndex/TwoElectronIntegralIndexTemplate.h"

#include "PTReferenceSpace.h"
#include "../../../app/CI/Selector/EnergyMap.h"

template <class MatrixType, class VectorType> class RepDiag;


template <class T> class Histogram;

template <class MatrixType>
class RepresentationMatrices : public CacheableObject {
public:
	
	RepresentationMatrices(TableKey key);
	~RepresentationMatrices();
	
	
	// copy constructor
	RepresentationMatrices(const RepresentationMatrices &);



//----------------------------------------------------------------------	

	
	INT	getP() const;
	INT	getNumberOfP3Mats() const;
	
	const RepresentationMatrix<MatrixType> *	getCbMat() const;
	const RepresentationMatrix<MatrixType> *	getExP3Mat() const;
	
//	RepresentationMatrix & operator [] (TableCase<MOType> &);

	// pointer to representation matrix to be mulplied with coulomb integrals
	const MatrixType *	getMatCoulombStart() const;	

	// pointer to representation matrix to be mulplied with exchange integrals
	const MatrixType *	getMatExchangeStart() const;
	
	// pointer to representation matrices for P=3 case
	const MatrixType *	getMatP3Start() const;
	
	// pointer to representation matrices for P=5 case
	const MatrixType *	getMatP5Start() const;
	
	
//----------------------------------------------------------------------	
	
	INT getOccupiedMemory() const;

//----------------------------------------------------------------------	

	INT	getNumberOfRows() const;
	INT	getNumberOfColumns() const;

//----------------------------------------------------------------------	
	
	void	printMats(ostream & s) const;
	

//----------------------------------------------------------------------	


static RepresentationMatrixFortranInterface<MatrixType> *repMatFInt;


protected:
	RepresentationMatrices();
	void	calcDense(MatrixType *buf);
	


INT	occupiedMemory;					// occupied memory of this class plus
									// matrix elements

typename RepresentationMatrix<MatrixType>::TMatDim 	matDim;


INT	P;								// P-case
INT	P3Mats;							// number of P=3 matrices

RepresentationMatrix<MatrixType>	*CbMat;		// pointer to matrix
RepresentationMatrix<MatrixType>	*ExP3Mat;	// pointer to matrix

};

template <class MatrixType> ostream & operator <<
		(ostream & s, const RepresentationMatrices<MatrixType> &);	



#include <math.h>


template <class MatrixType>
inline
INT RepresentationMatrices<MatrixType>::getP() const
{	return P;	}

template <class MatrixType>
inline
INT RepresentationMatrices<MatrixType>::getNumberOfP3Mats() const
{	return P3Mats;	}

template <class MatrixType>
inline
const RepresentationMatrix<MatrixType> *	RepresentationMatrices<MatrixType>::getCbMat() const
{	return	CbMat;	}

template <class MatrixType>
inline
const RepresentationMatrix<MatrixType> *	RepresentationMatrices<MatrixType>::getExP3Mat() const
{	return	ExP3Mat;	}
	
template <class MatrixType>
inline
INT RepresentationMatrices<MatrixType>::getOccupiedMemory() const
{	return occupiedMemory;	}

template <class MatrixType>
inline
const MatrixType *	RepresentationMatrices<MatrixType>::getMatCoulombStart() const
{	return CbMat->getP();	}

template <class MatrixType>
inline
const MatrixType *	RepresentationMatrices<MatrixType>::getMatExchangeStart() const
{	return ExP3Mat->getP();	}

template <class MatrixType>
inline
const MatrixType *	RepresentationMatrices<MatrixType>::getMatP3Start() const
{	return ExP3Mat->getP();	}

template <class MatrixType>
inline
const MatrixType *	RepresentationMatrices<MatrixType>::getMatP5Start() const
{	return ExP3Mat->getP();	}


template <class MatrixType>
inline
INT	RepresentationMatrices<MatrixType>::getNumberOfRows() const
{	return matDim.rows;	}

template <class MatrixType>
inline
INT	RepresentationMatrices<MatrixType>::getNumberOfColumns() const
{	return matDim.cols;	}

template <class MatrixType>
ostream &	operator<<(ostream &s, const RepresentationMatrices<MatrixType> & repMats)
{
	s << "P=" << repMats.getP() << endl << endl;
	repMats.printMats(s);
	return s;
}


#endif
