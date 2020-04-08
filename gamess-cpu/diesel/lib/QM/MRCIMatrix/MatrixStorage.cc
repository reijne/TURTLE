//***********************************************************************
//
//	Name:			MatrixStorage.cc
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

#include "MatrixStorage.h"

#include <iostream>
#include <iomanip>
#include <string>
#include <unistd.h>

using namespace std;

INT writeHamiltonOnly = 0;
template <class MatrixType, class VectorType> 
const INT MatrixStorage<MatrixType, VectorType>::blocksizeInts;
template <class MatrixType, class VectorType> 
const INT MatrixStorage<MatrixType, VectorType>::blocksizeBytes;



template <class MatrixType, class VectorType>
MatrixStorage<MatrixType, VectorType>::MatrixStorage()
{
	dim = 0;
	nEntries = 0;
	p = NULL;
	fio = NULL;

	wo = 0;

	if ( writeHamiltonOnly )
	{
		fio = new FortranFileIO("MatrixStorageWrite.dat");

		fio->write(&dim, sizeof(dim));
		fio->write((INT *)&blocksizeInts, sizeof(blocksizeInts));
	}

}


template <class MatrixType, class VectorType>
MatrixStorage<MatrixType, VectorType>::MatrixStorage(INT _dim)
{
	dim = _dim;
	nEntries = 0;
	p = NULL;
	fio = NULL;


	wo = 0;

	if ( writeHamiltonOnly )
	{
		unlink("MatrixStorageWrite.dat");
		fio = new FortranFileIO("MatrixStorageWrite.dat");

		fio->write(&dim, sizeof(dim));
		fio->write((INT *) &blocksizeInts, sizeof(blocksizeInts));
	}
}


template <class MatrixType, class VectorType>
MatrixStorage<MatrixType, VectorType>::
	MatrixStorage(const MatrixStorage<MatrixType, VectorType> & mat)
{
	dim = mat.dim;
	nEntries = mat.nEntries;
	p = new TEntry[nEntries];
	memcpy(p, mat.p, nEntries*sizeof(TEntry));
}


template <class MatrixType, class VectorType>
MatrixStorage<MatrixType, VectorType> & MatrixStorage<MatrixType, VectorType>::operator = (
	const MatrixStorage<MatrixType, VectorType> & mat)
{
	if ( p )
		delete p;

	dim = mat.dim;
	nEntries = mat.nEntries;
	p = new TEntry[nEntries];
	memcpy(p, mat.p, nEntries*sizeof(TEntry));
	
	return *this;
}


template <class MatrixType, class VectorType>
MatrixStorage<MatrixType, VectorType>::~MatrixStorage()
{
	if ( p )
		delete p;


	if ( fio )
	{
		memset(&buffer[wo], 0,
			(blocksizeInts-wo)*sizeof(TEntry));
		fio->write(buffer, blocksizeBytes);

		delete fio;
	}

}



template <class MatrixType, class VectorType>
void	MatrixStorage<MatrixType, VectorType>::mult(const VectorType *x, VectorType *y) const
{
	memset(y, 0, dim*sizeof(VectorType));

TEntry	*pp = p;	
	for ( INT i=0 ; i<nEntries ; i++ )
	{
		y[pp->row] += pp->v * x[pp->col];
		if ( pp->row != pp->col )
			y[pp->col] += pp->v * x[pp->row];
		pp++;
	}
}

template <class MatrixType, class VectorType>
void	MatrixStorage<MatrixType, VectorType>::mult(const VectorType *x, VectorType *y, INT n) const
{


VectorType* x1 = new VectorType[n*dim];
VectorType* y1 = new VectorType[n*dim];
VectorType	*xx1 = x1;
VectorType	*yy1 = y1;
	memset(y1, 0, n*dim*sizeof(VectorType));

	// reorganize array for faster access
	for ( INT i=0 ; i<dim ; i++ )
		for ( INT j=0 ; j<n ; j++ )
			*xx1++ = x[j*dim+i];


TEntry	*pp = p;	
	for ( INT i=0 ; i<nEntries ; i++ )
	{
	INT	f = (pp->row != pp->col);

	VectorType	*xxr = x1 + pp->row*n;
	VectorType	*xxc = x1 + pp->col*n;
	VectorType	*yyr = y1 + pp->row*n;
	VectorType	*yyc = y1 + pp->col*n;
	VectorType v = pp->v;
		
		for ( INT j=0 ; j<n ; j++ )
		{
			*yyr++ +=     v * *xxc++;
			*yyc++ += f * v * *xxr++;
		}
		pp++;
	}


	// reorder array for succeeding processing
	yy1 = y1;
	for ( INT i=0 ; i<dim ; i++ )
		for ( INT j=0 ; j<n ; j++ )
			y[j*dim+i] = *yy1++;
	delete[] y1;
	delete[] x1;
}



template <class MatrixType, class VectorType>
double	MatrixStorage<MatrixType, VectorType>::getSparsity() const
{
	return (1.0 - 1.0*nEntries/(dim*(dim+1)/2));
}

template <class MatrixType, class VectorType>
LONG_LONG_INT MatrixStorage<MatrixType, VectorType>::getMemory(LONG_LONG_INT n)
{
	return n*sizeof(TEntry);
}


template <class MatrixType, class VectorType>
ostream & operator << (ostream &os, const MatrixStorage<MatrixType, VectorType> &m)
{
MatrixType	*row = new MatrixType[m.dim];

	for ( INT i=0 ; i<m.dim ; i++ )
	{
		memset(row, 0, i*sizeof(MatrixType));
		
		for ( INT j=0 ; j<m.nEntries ; j++ )
			if ( m.p[j].row==i )
				row[m.p[j].col] = m.p[j].v;

		for ( INT j=0 ; j<=i ; j++ )
			os << setw(12) << row[j];
		os << endl;
	}
	
	
	delete row;

	return os;
}



template class MatrixStorage<float, float>;
template class MatrixStorage<float, double>;
template class MatrixStorage<double, float>;
template class MatrixStorage<double, double>;





