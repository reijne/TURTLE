//***********************************************************************
//
//	Name:			SLEMatrix.h
//
//	Description:	base class for matrices in solution of linear equations
//					solve A*x=b for x, A: VERY large matrix
//
//	Author:			Michael Hanrath
//
//	Version:		0.0
//
//	Date:			16.09.1998
//
//
//
//
//***********************************************************************


#ifndef __SLEMatrix_h
#define __SLEMatrix_h

#include "../../../config.h"


#include "IterationMatrix.h"


template <class VectorType> class BufferedVector;

template <class MatrixType, class VectorType>
class SLEMatrix : 
	public IterationMatrix<MatrixType, VectorType> {
public:
	SLEMatrix(INT totalDim) : IterationMatrix<MatrixType, VectorType>(totalDim) {};
	virtual ~SLEMatrix() {};



	

	// calculate resiual vector
	//          x  = x - b
	virtual void	calcResidual(BufferedVector<VectorType> &x) const = 0;


	// perform multiplication with invers of diagonal
	//                  1
	//          x  = -------  x 
	//           i      A      i
	//                   ii
	virtual void	multInvDiag(BufferedVector<VectorType> &x) const = 0;


	// perform multiplication on n vectors
	// y = A*x
	virtual void	mult(const BufferedVector<VectorType> &x, BufferedVector<VectorType> &y) const = 0;


protected:
	SLEMatrix() : IterationMatrix<MatrixType, VectorType>(0) {};
};





#endif
