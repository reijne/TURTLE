//***********************************************************************
//
//	Name:			DavidsonMatrixStorage.cc
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

#include "DavidsonMatrixStorage.h"

#include <iostream>
#include <iomanip>

#include "../MRTree/Diag/NExternalsDiag.h"

using namespace std;

template <class MatrixType, class VectorType>
DavidsonMatrixStorage<MatrixType, VectorType>::
	DavidsonMatrixStorage(const MRCIMatrix<MatrixType, VectorType> *_mrciMatrix) :
	MatrixStorage<MatrixType, VectorType>(_mrciMatrix->getNExternalsDiag()->
		getNumberOfTotalSpinAdaptedFunctions()),
	DavidsonMatrix<MatrixType, VectorType>(_mrciMatrix->getNExternalsDiag()->
		getNumberOfRefConfSpinAdaptedFunctions(),
		_mrciMatrix->getNExternalsDiag()->
		getNumberOfTotalSpinAdaptedFunctions())
{
	mrciMatrix = _mrciMatrix;


	cout << "storing Hamilton Matrix..." << endl;

	refMatDim = mrciMatrix->getNExternalsDiag()->
		getNumberOfRefConfSpinAdaptedFunctions();
	referenceMatrix = new MatrixType[refMatDim*refMatDim];
	
	mrciMatrix->getRefMatrix(referenceMatrix);
	
	mrciMatrix->getTotalMatrix(*this);
	
//	cout << *this << endl;
	
	cout << "sparsity = " << this->getSparsity()*100 << "%" << endl;
}


template <class MatrixType, class VectorType>
DavidsonMatrixStorage<MatrixType, VectorType>::~DavidsonMatrixStorage()
{
	delete referenceMatrix;
}



template <class MatrixType, class VectorType>
void	DavidsonMatrixStorage<MatrixType, VectorType>::getDiagonal(MatrixType *pp) const
{
	memset(pp, 0, this->totalDim*sizeof(MatrixType));

        typename DavidsonMatrixStorage<MatrixType,VectorType>::TEntry *pe = this->p;

	for ( INT i=0 ; i<this->nEntries ; i++ )
	{
		if ( pe->col==pe->row )
			pp[pe->col] = pe->v;
		pe++;
	}
}

template <class MatrixType, class VectorType>
void	DavidsonMatrixStorage<MatrixType, VectorType>::getRefMatrix(MatrixType *pp) const
{
	memcpy(pp, referenceMatrix, refMatDim*refMatDim*sizeof(MatrixType));
}






template <class MatrixType, class VectorType>
void	DavidsonMatrixStorage<MatrixType, VectorType>::getTotalMatrix(
	MatrixType *pp) const
{

	for ( INT k=0 ; k<this->nEntries ; k++ )
	{
	INT	i = this->p[k].row;
		i = i*(i+1)/2;
		pp[i + this->p[k].col] = this->p[k].v;
	}
}


template <class MatrixType, class VectorType>
void	DavidsonMatrixStorage<MatrixType, VectorType>::getTotalMatrix(MatrixStorage<MatrixType, VectorType> &m) const
{
}

template <class MatrixType, class VectorType>
void	DavidsonMatrixStorage<MatrixType, VectorType>::
	mult(const VectorType *x, VectorType *y, INT n) const
{
/*	for ( INT k=0 ; k<n ; k++ )
		MatrixStorage<MatrixType, VectorType>::mult(x+k*totalDim, y+k*totalDim);
	return;
*/	
		if ( n==1 )
		MatrixStorage<MatrixType, VectorType>::mult(x, y);
	else
		MatrixStorage<MatrixType, VectorType>::mult(x, y, n);
}


template <class MatrixType, class VectorType>
ostream & operator<< (ostream &os, const DavidsonMatrixStorage<MatrixType, VectorType> &m)
{
MatrixType	*row = new MatrixType[m.dim];

	for ( INT i=0 ; i<m.dim ; i++ )
	{
		memset(row, 0, m.dim*sizeof(MatrixType));
		
		for ( INT j=0 ; j<m.nEntries ; j++ )
		{
			if ( m.p[j].row==i )
				row[m.p[j].col] = m.p[j].v;
			if ( m.p[j].col==i )
				row[m.p[j].row] = m.p[j].v;
		}

		for ( INT j=0 ; j<m.dim ; j++ )
			os << setprecision(8) << setw(18) << row[j];
		os << endl;
	}
	
	
	delete row;

	return os;
}


template class DavidsonMatrixStorage<float, float>;
template class DavidsonMatrixStorage<float, double>;
template class DavidsonMatrixStorage<double, float>;
template class DavidsonMatrixStorage<double, double>;
