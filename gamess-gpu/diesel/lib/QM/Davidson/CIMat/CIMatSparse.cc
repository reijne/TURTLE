//***********************************************************************
//
//	Name:			CIMatSparse.cc
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

#include "CIMatSparse.h"
#include "MatSel.h"

#include "../../../Container/DiskBuffer.h"

#include "../DavidsonMatrix.h"

#include <stdlib.h>
#include <math.h>
#include <string>
#include <iostream>
#include <iomanip>

using namespace std;

CIMatSparse::CIMatSparse(INT refDim, INT totalDim) : CIMat(refDim, totalDim)
{
}

CIMatSparse::CIMatSparse(istream &is) : CIMat(0, 0)
{
	is >> refDim;
	refMat = new MatrixType[refDim*(refDim+1)/2];
	memset(refMat, 0, refDim*(refDim+1)/2*sizeof(MatrixType));
	
	for ( INT i=0 ; i<refDim*(refDim+1)/2 ; i++ )
		is >> refMat[i];



	is >> totalDim;
	
	for ( INT i=0 ; i<totalDim ; i++ )
		for ( INT j=0 ; j<=i ; j++ )
		{
		double	v;
			is >> v;
			addEntry(i, j, v);
		}
}




CIMatSparse::~CIMatSparse()
{
}


CIMatSparse::CIMatSparse(const CIMatSparse &mat) :
	CIMat(mat), MatrixStorage<MatrixType, VectorType>(mat)
{
}

CIMatSparse::CIMatSparse(const CIMatSparse &mat, INT _totalDim) :
	CIMat(mat)
{
	totalDim = _totalDim;
	for ( INT i=0 ; i<nEntries ; i++ )
		if ( mat.p[i].row<totalDim && mat.p[i].col<totalDim )
			addEntry(mat.p[i].row, mat.p[i].col, mat.p[i].v);
}


CIMatSparse::CIMatSparse(const CIMatSparse &mat, const MatSel &sel) : 
	CIMat(mat.refDim, sel.getDim())

{	
	for ( INT i=0 ; i<totalDim ; i++ )
		for ( INT j=i ; j<totalDim ; j++ )
			(*this)(i, j) = mat(sel[i], sel[j]);
}


CIMatSparse & CIMatSparse::operator = (const CIMatSparse & mat)
{
	(MatrixStorage<MatrixType, VectorType> &) *this = mat;
	(CIMat &) *this = mat;

	return *this;
}

MatrixType CIMatSparse::zero = 0;	


MatrixType &	CIMatSparse::operator()(INT i, INT j)
{
	for ( INT k=0 ; k<nEntries ; k++ )
		if ( p[k].row==i && p[k].col==j  ||  p[k].row==j && p[k].col==i )
			return p[k].v;

//	addEntry(i, j, 0);
	
//	return p[nEntries-1].v;
	return zero;
}

MatrixType	CIMatSparse::operator()(INT i, INT j) const
{
	for ( INT k=0 ; k<nEntries ; k++ )
		if ( p[k].row==i && p[k].col==j  ||  p[k].row==j && p[k].col==i )
			return p[k].v;

	return 0;
}


void	CIMatSparse::getTotalMatrix(MatrixType *pp) const
{

	for ( INT k=0 ; k<nEntries ; k++ )
	{
	INT	i = p[k].row;
		i = i*(i+1)/2;
		pp[i + p[k].col] = p[k].v;
	}
}


void	CIMatSparse::getTotalMatrix(MatrixStorage<MatrixType, VectorType> &m) const
{
}

void	CIMatSparse::getDiagonal(MatrixType *pp) const
{
	for ( INT k=0 ; k<nEntries ; k++ )
	{
		if ( p[k].row == p[k].col )
			pp[p[k].row] = p[k].v;
	}
}


void	CIMatSparse::mult(const double *x, double *y, INT n) const
{
	for ( INT k=0 ; k<n ; k++ )
		MatrixStorage<MatrixType, VectorType>::mult(x+k*totalDim, y+k*totalDim);
}


