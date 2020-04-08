//***********************************************************************
//
//	Name:			DavidsonMatrix.cc
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

#include "DavidsonMatrix.h"
#include "../../Container/DiskBuffer.h"

#include <string>
#include <iostream>

using namespace std;

template <class MatrixType, class VectorType>
void	DavidsonMatrix<MatrixType, VectorType>::
	mult(DiskBuffer *xBuf, DiskBuffer *yBuf, INT start, INT end) const
{

INT	n = end - start + 1;
VectorType	*x = new VectorType[this->totalDim*n];
VectorType	*y = new VectorType[this->totalDim*n];


	for ( INT j=0 ; j<n ; j++ )
		xBuf->get(start + j, x+j*this->totalDim);

	memset(y, 0, this->totalDim*n*sizeof(MatrixType));

	mult(x, y, n);

/*
	cout << endl;
	for ( INT i=0 ; i<this->totalDim ; i++ )
		cout << x[i] << "\t";
	cout << endl;
	for ( INT i=0 ; i<this->totalDim ; i++ )
		cout << y[i] << "\t";
	cout << endl;
	cout << endl;
*/

	for ( INT j=0 ; j<n ; j++ )
		yBuf->put(start + j, y+j*this->totalDim);

	delete x;
	delete y;
}

template class DavidsonMatrix<float, float>;
template class DavidsonMatrix<float, double>;
template class DavidsonMatrix<double, float>;
template class DavidsonMatrix<double, double>;
