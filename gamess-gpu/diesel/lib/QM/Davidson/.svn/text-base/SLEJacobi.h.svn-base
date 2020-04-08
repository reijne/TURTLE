//***********************************************************************
//
//	Name:			SLEJacobi.h
//
//	Description:	solution of large linear equation system
//					by Jacobi's method 
//					(total step method: 
//						inverse of A is approximated
//						by inverse of diagonal)
//
//	Author:			Michael Hanrath
//
//	Version:		0.0
//
//	Date:			16.08.1998
//
//
//
//
//***********************************************************************


#ifndef __SLEJacobi_h
#define __SLEJacobi_h

#include "../../../config.h"


template <class MatrixType, class VectorType> class SLEMatrix;

template <class VectorType> class BufferedVector;


template <class MatrixType, class VectorType>
class SLEJacobi {
public:
	SLEJacobi(const SLEMatrix<MatrixType, VectorType> *);
	~SLEJacobi();
	

	void	iterate();

	const	BufferedVector<VectorType>	&getX() const;
	const	BufferedVector<VectorType>	&getB() const;


private:
const SLEMatrix<MatrixType, VectorType> *sleMatrix;
BufferedVector<VectorType>	&x;
BufferedVector<VectorType>	&b;

static const VectorType	epsilon;
};





#endif
