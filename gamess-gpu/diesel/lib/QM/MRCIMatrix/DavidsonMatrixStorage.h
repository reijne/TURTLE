//***********************************************************************
//
//	Name:			DavidsonMatrixStorage.h
//
//	Description:	
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
//
//***********************************************************************

#ifndef __DavidsonMatrixStorage_h
#define __DavidsonMatrixStorage_h

#include <iostream>
using std::ostream;

#include "../../../config.h"

#include "MatrixStorage.h"
#include "../Davidson/DavidsonMatrix.h"

template <class MatrixType, class VectorType> class MRCIMatrix;

template <class MatrixType, class VectorType> class DavidsonMatrixStorage;

template <class MatrixType, class VectorType> 
       ostream & operator<< (ostream &, const DavidsonMatrixStorage<MatrixType, VectorType> &);

template <class MatrixType, class VectorType>
class DavidsonMatrixStorage : 
	public MatrixStorage<MatrixType, VectorType>,
	public DavidsonMatrix<MatrixType, VectorType> {

public:
	DavidsonMatrixStorage(const MRCIMatrix<MatrixType, VectorType> *);
	~DavidsonMatrixStorage();
	

//----------------------------------------------------------------------------
//	functions from DavidsonMatrix

	// get diagonal of matrix
	void	getDiagonal(MatrixType *) const;
	
	// get reference part of matrix
	void	getRefMatrix(MatrixType *) const;

	// get total matrix (attention: probably VERY large)
	void	getTotalMatrix(MatrixType *) const;
	void	getTotalMatrix(MatrixStorage<MatrixType, VectorType> &) const;

	// get index of nth reference in total matrix
	INT	getRefIndex(INT) const;

	// perform multiplication 
	// y = A*x
	void	mult(const VectorType *x, VectorType *y, INT n) const;


//----------------------------------------------------------------------------

	friend ostream & operator<< <> (ostream &, const DavidsonMatrixStorage &);

private:
const MRCIMatrix<MatrixType, VectorType>	*mrciMatrix;
INT	refMatDim;
MatrixType	*referenceMatrix;
};


#include "MRCIMatrix.h"


template <class MatrixType, class VectorType>
inline
INT	DavidsonMatrixStorage<MatrixType, VectorType>::getRefIndex(INT i) const
{	return mrciMatrix->getRefIndex(i);	}


#endif
