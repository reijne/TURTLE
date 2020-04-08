//***********************************************************************
//
//	Name:			ConjGrad.h
//
//	Description:	solution of large linear equation system
//					by method of conjugated gradient
//
//	Author:			Michael Hanrath
//
//	Version:		0.0
//
//	Date:			30.08.1998
//
//	Ref:			
//
//
//***********************************************************************


#ifndef __ConjGrad_h
#define __ConjGrad_h

#include "../../../config.h"

template <class MatrixType, class VectorType> class SLEMatrix;

template <class VectorType> class BufferedVector;


template <class MatrixType, class VectorType>
class ConjGrad {
public:
	ConjGrad(const SLEMatrix<MatrixType, VectorType> *);
	~ConjGrad();
	

	void	iterate();
	void	iterateS();

	const	BufferedVector<VectorType>	& getX() const;
	const	BufferedVector<VectorType>	& getB() const;
	
private:
const SLEMatrix<MatrixType, VectorType> *sleMatrix;
BufferedVector<VectorType>	&x;
BufferedVector<VectorType>	&b;

static const VectorType	epsilon;
};





#endif
