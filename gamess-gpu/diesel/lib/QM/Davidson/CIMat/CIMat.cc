//***********************************************************************
//
//	Name:			CIMat.cc
//
//	Description:	
//
//	Author:			Michael Hanrath
//
//	Version:		0.0
//
//	Date:			17. Dec 1998
//
//***********************************************************************

#include "CIMat.h"
#include "MatSel.h"

#include "../../../Container/DiskBuffer.h"

#include <stdlib.h>
#include <math.h>
#include <string>
#include <iostream>
#include <iomanip>

using namespace std;

CIMat::CIMat(istream &is) :
		DavidsonMatrix<MatrixType, VectorType>(0, 0)
{
	is >> refDim;
	refMat = new MatrixType[refDim*refDim];
	memset(refMat, 0, refDim*refDim*sizeof(MatrixType));
	
	for ( INT i=0 ; i<refDim*refDim ; i++ )
		is >> refMat[i];

	index = new INT[refDim];
	memset(index, 0, refDim*sizeof(INT));
	for ( INT i=0 ; i<refDim ; i++ )
		is >> index[i];

}

CIMat::CIMat(const CIMat &mat) :
	DavidsonMatrix<MatrixType, VectorType>(mat.getRefDim(), mat.getTotalDim())
{
	refMat = new MatrixType[refDim*refDim];
	memcpy(refMat, mat.refMat, refDim*refDim*sizeof(MatrixType));

	index = new INT[refDim];
	memcpy(index, mat.index, refDim*sizeof(INT));
}

CIMat::CIMat(const CIMat &mat, INT totalDim) : 
		DavidsonMatrix<MatrixType, VectorType>(mat.getRefDim(), totalDim)
{
	refMat = new MatrixType[refDim*refDim];
	memcpy(refMat, mat.refMat, refDim*refDim*sizeof(MatrixType));

	index = new INT[refDim];
	memcpy(index, mat.index, refDim*sizeof(INT));
}


CIMat & CIMat::operator = (const CIMat & mat)
{
	if ( refMat )
		delete refMat;
		
	refMat = new MatrixType[refDim*refDim];
	memcpy(refMat, mat.refMat, refDim*refDim*sizeof(MatrixType));


	if ( index )
		delete index;
		
	index = new INT[refDim];
	memcpy(index, mat.index, refDim*sizeof(INT));

	return *this;
}


CIMat::CIMat(INT refDim, INT dim) :
		DavidsonMatrix<MatrixType, VectorType>(refDim, dim)
{
	refMat = new MatrixType[refDim*refDim];
	memset(refMat, 0, refDim*refDim*sizeof(MatrixType));

	index = new INT[refDim];
	memset(index, 0, refDim*sizeof(INT));
}


CIMat::~CIMat()
{
	if ( refMat )
		delete refMat;
		
	if ( index )
		delete index;
}

void	CIMat::getRefMatrix(MatrixType *pp) const
{
	memcpy(pp, refMat, refDim*refDim*sizeof(MatrixType));
}



INT	CIMat::getRefIndex(INT i) const
{	return index[i];	}


void	CIMat::setRandomCI(double dens,
		double diag, double diagDispl, double nonDiagDispl, INT randSeed)
{
	srand(randSeed);

	for ( INT i=0 ; i<totalDim ; i++ )
		for ( INT j=i ; j<totalDim ; j++ )
		{
			if ( 1.0*rand() / RAND_MAX >dens )
				(*this)(i,j) = nonDiagDispl*(1.0*rand() / RAND_MAX - 0.5);

			if ( i==j )
				(*this)(i,j) = diag + diagDispl*(1.0*rand() / RAND_MAX - 0.5);
			
		}
}


void	CIMat::getDiagonal(MatrixType *pp) const
{
	for ( INT i=0 ; i<totalDim ; i++ )
		pp[i] = (*this)(i, i);
}

struct TDat {
	MatrixType	v;
	INT	i;
	};

static int cmp(const void *p1, const void *p2)
{
	if ( static_cast<const TDat *>(p1)->v>static_cast<const TDat *>(p2)->v )
		return 1;
	if ( static_cast<const TDat *>(p1)->v<static_cast<const TDat *>(p2)->v )
		return -1;
	return 0;
}

void	CIMat::chooseRefMat()
{
	memset(refMat, 0, refDim*refDim*sizeof(MatrixType));
	
TDat* d = new TDat[totalDim];
	
/*	cout << ":::::::::::::::::::::::::::" << endl;
	cout << *this << endl;
	cout << ":::::::::::::::::::::::::::" << endl;
*/
MatrixType* v = new MatrixType[totalDim];

	getDiagonal(v);

	for ( INT i=0 ; i<totalDim ; i++ )
	{
		d[i].v = v[i];
		d[i].i = i;
	}
	qsort(d, totalDim, sizeof(TDat), cmp);


	for ( INT i=0 ; i<refDim ; i++ )
		index[i] = d[i].i;


	for ( INT i=0 ; i<refDim ; i++ )
		for ( INT j=0 ; j<refDim ; j++ )
		{
//			cout << i << " " << j << " " << index[i] << " " << index[j] << endl;
			refMat[refDim*i + j] = (*this)(index[i], index[j]);
		}
        delete v;
	delete d;
}


void	CIMat::mult(const double *x, double *y, INT n) const
{
	for ( INT i=0 ; i<totalDim ; i++ )
		for ( INT j=0 ; j<totalDim ; j++ )
			for ( INT k=0 ; k<n ; k++ )
				y[k*totalDim + i] += (*this)(i, j)*x[k*totalDim + j];
}






ostream & operator << (ostream & os, const CIMat &mat) 
{
	for ( INT i=0 ; i<mat.getTotalDim() ; i++ )
	{
		for ( INT j=0 ; j<mat.getTotalDim() ; j++ )
			os << setw(16) << ((CIMat &) mat)(i, j);
		
		os << endl;
	}
	return os;
}
