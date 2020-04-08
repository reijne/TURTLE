//***********************************************************************
//
//	Name:			IterationMatrix.h
//
//	Description:	base class for matrices used in
//					- solution of linear equations
//					- davidson iteration
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


#ifndef __IterationMatrix_h
#define __IterationMatrix_h


#include "../../../config.h"


template <class MatrixType, class VectorType> class MatrixStorage;

template <class MatrixType, class VectorType>
class IterationMatrix {
public:
	IterationMatrix(INT totalDim) : totalDim(totalDim) {};
	virtual ~IterationMatrix() {};

	INT getTotalDim() const;



	// get diagonal of matrix
	virtual void	getDiagonal(MatrixType *) const = 0;


	// get total matrix (attention: probably VERY large)
	virtual void	getTotalMatrix(MatrixType *) const = 0;
	virtual void	getTotalMatrix(
		MatrixStorage<MatrixType, VectorType> &) const = 0;

protected:
	
INT	totalDim;			// total dimension
};


template <class MatrixType, class VectorType>
inline
INT IterationMatrix<MatrixType, VectorType>::getTotalDim() const
{	return totalDim;	}



#endif
