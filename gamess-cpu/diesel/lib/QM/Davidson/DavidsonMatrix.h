//***********************************************************************
//
//	Name:			DavidsonMatrix.h
//
//	Description:	base class for matrices used in davidson iteration
//
//	Author:			Michael Hanrath
//
//	Version:		0.0
//
//	Date:			14.08.1998
//
//
//
//
//***********************************************************************


#ifndef __DavidsonMatrix_h
#define __DavidsonMatrix_h

#include "../../../config.h"


#include "IterationMatrix.h"



class DiskBuffer;

template <class MatrixType, class VectorType>
class DavidsonMatrix : 
	public IterationMatrix<MatrixType, VectorType> {
public:
	DavidsonMatrix(INT refDim, INT totalDim) :
		IterationMatrix<MatrixType, VectorType>(totalDim), refDim(refDim) {};
	virtual ~DavidsonMatrix() {};

	INT getRefDim() const;

	
	// get reference part of matrix
	virtual void	getRefMatrix(MatrixType *) const = 0;


	// get index of nth reference in total matrix
	virtual INT	getRefIndex(INT) const = 0;

	// perform multiplication 
	// y = A*x
	virtual void	mult(DiskBuffer *x, DiskBuffer *y, INT start, INT end) const;

protected:
	virtual void	mult(const VectorType *x, VectorType *y, INT n) const = 0;
	
INT	refDim;				// dimension of reference part

};


template <class MatrixType, class VectorType>
inline
INT DavidsonMatrix<MatrixType, VectorType>::getRefDim() const
{	return refDim;	}



#endif
