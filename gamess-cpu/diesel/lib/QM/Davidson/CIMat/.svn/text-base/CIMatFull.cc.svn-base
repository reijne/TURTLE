//***********************************************************************
//
//	Name:			CIMatFull.cc
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

#include "CIMatFull.h"
#include "MatSel.h"

#include "../../../Container/DiskBuffer.h"

#include <stdlib.h>
#include <math.h>
#include <string>
#include <iostream>
#include <iomanip>

using namespace std;

CIMatFull::CIMatFull(INT refDim, INT totalDim) :
	CIMat(refDim, totalDim)
{
	p = new MatrixType[totalDim*(totalDim+1)/2];
	memset(p, 0, totalDim*(totalDim+1)/2*sizeof(MatrixType));
}


CIMatFull::CIMatFull(istream &is) : CIMat(is)
{

	is >> totalDim;
	p = new MatrixType[totalDim*(totalDim+1)/2];
	memset(p, 0, totalDim*(totalDim+1)/2*sizeof(MatrixType));
	
	for ( INT i=0 ; i<totalDim*(totalDim+1)/2 ; i++ )
		is >> p[i];
}


CIMatFull::CIMatFull(const MatrixType *_p, INT refDim, INT dim) :
	CIMat(refDim, dim)
{
	p = new MatrixType[totalDim*(totalDim+1)/2];
	memcpy(p, _p, totalDim*(totalDim+1)/2*sizeof(MatrixType));
	chooseRefMat();
}



CIMatFull::~CIMatFull()
{
	if ( p )
		delete p;
}


CIMatFull::CIMatFull(const CIMatFull &mat) :
	CIMat(mat)
{
	p = new MatrixType[totalDim*(totalDim+1)/2];
	memcpy(p, mat.p, totalDim*(totalDim+1)/2*sizeof(MatrixType));
}

CIMatFull::CIMatFull(const CIMatFull &mat, INT totalDim) : 
	CIMat(mat, totalDim)
{
	p = new MatrixType[totalDim*(totalDim+1)/2];
	memcpy(p, mat.p, totalDim*(totalDim+1)/2*sizeof(MatrixType));
}


CIMatFull::CIMatFull(const CIMatFull &mat, const MatSel &sel) : 
	CIMat(mat.refDim, sel.getDim())
{
	refMat = new MatrixType[refDim*refDim];
	memcpy(refMat, mat.refMat, refDim*refDim*sizeof(MatrixType));
	
	
	
	index = new INT[refDim];
	for ( INT i=0 ; i<refDim ; i++ )
	{
	INT	j;
		for ( j=0 ; j<sel.getDim() ; j++ )
			if ( sel[j]==mat.index[i] )
			{
//				cout << i << " " << j << " " << mat.index[i] << " " << sel[j] << endl;
				index[i] = j;
				break;
			}
//			index[i] = mat.index[i];
	}

	p = new MatrixType[totalDim*(totalDim+1)/2];
	
	for ( INT i=0 ; i<totalDim ; i++ )
		for ( INT j=i ; j<totalDim ; j++ )
			(*this)(i, j) = mat(sel[i], sel[j]);
}


CIMatFull & CIMatFull::operator = (const CIMatFull & mat)
{
	static_cast<CIMat &>(*this) = mat;
	if ( p )
		delete p;
		
	totalDim = mat.totalDim;
	p = new MatrixType[totalDim*(totalDim+1)/2];
	memcpy(p, mat.p, totalDim*(totalDim+1)/2*sizeof(MatrixType));
	return *this;
}


void	CIMatFull::getTotalMatrix(MatrixType *pp) const
{
	memcpy(pp, p, totalDim*(totalDim+1)/2*sizeof(MatrixType));
}

void	CIMatFull::getTotalMatrix(MatrixStorage<MatrixType, VectorType> &) const
{
}







