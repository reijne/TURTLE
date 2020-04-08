//***********************************************************************
//
//	Name:			CIEig.cc
//
//	Description:	
//
//	Author:			Michael Hanrath
//
//	Version:		0.0
//
//	Date:			15. Aug 1998
//
//***********************************************************************

#include "CIEig.h"
#include "CIMat.h"


#include "../Davidson.h"
#include "../Roots.h"

#include <math.h>
#include <string>
#include <iostream>
#include <iomanip>

#include "../../../Math/FortranNumeric/FortranEigenProblems.h"

using namespace std;

CIEig::CIEig(const CIMat &mat) :
	nRoots(0)
{
	dim = mat.getTotalDim();

INT	nV = dim;
	pE = new double[dim];
	memset(pE, 0, dim*sizeof(double));
	pV = new double[dim*dim];
	memset(pV, 0, dim*dim*sizeof(double));

double	*pMat = new double[dim*dim];
	mat.getTotalMatrix(pMat);

	hqrii1(
		&dim, &dim, pMat, pE, pV, &nV
	);

	for ( INT i=0 ; i<nV ; i++ )
	{
	double	*p = pV+i*dim;
	double	sum = 0;
		for ( INT j=0 ; j<dim ; j++ )
			sum += *p * *p++;
		sum = 1/sqrt(sum);
		p = pV+i*dim;
		for ( INT j=0 ; j<dim ; j++ )
			*p++ *= sum;
	}
	delete pMat;
}



CIEig::CIEig(const CIMat &mat, INT nRoots) :
	nRoots(nRoots)
{
	cout << "nRoots=" << nRoots << endl;
	dim = mat.getTotalDim();


	pE = new double[nRoots];
	memset(pE, 0, nRoots*sizeof(double));
	pV = new double[nRoots*dim];
	memset(pV, 0, nRoots*dim*sizeof(double));

Roots	roots(nRoots);
	cout << roots;

Davidson<MatrixType, VectorType>	davidson(dim, mat.getRefDim(), roots, mat, 0, 0);
	davidson.start();


	davidson.iterate(
		100,
		Davidson<MatrixType, VectorType>::Energy,
		1e-8);


	for ( INT i=0 ; i<nRoots ; i++ )
		pE[i] = roots.getRoot(i);

DiskBuffer	*evBuf = new DiskBuffer(dim*sizeof(VectorType),
	"Eigenvectors.dat", DiskBuffer::noTempDir);
	for ( INT i=0 ; i<nRoots ; i++ )
		evBuf->get(i, pV + i*dim);


	delete evBuf;
/*			for ( INT i=0 ; i<nRoots ; i++ )
	{
	double	*p = pV+i*dim;
	double	sum = 0;
		for ( INT j=0 ; j<dim ; j++ )
			sum += *p * *p++;
		sum = 1/sqrt(sum);
		p = pV+i*dim;
		for ( INT j=0 ; j<dim ; j++ )
			*p++ *= sum;
	}
*/
}



CIEig::~CIEig()
{
	if ( pE )
		delete pE;
	if ( pV )
		delete pV;
}

CIEig::CIEig(const CIEig &eig)
{
	dim = eig.getDim();

	pE = new double[dim];
	memcpy(pE, eig.pE, dim*sizeof(double));
	pV = new double[dim*dim];
	memcpy(pV, eig.pV, dim*dim*sizeof(double));
}


CIEig & CIEig::operator = (const CIEig &eig)
{
	if ( pE )
		delete pE;
	if ( pV )
		delete pV;

	dim = eig.getDim();

	pE = new double[dim];
	memcpy(pE, eig.pE, dim*sizeof(double));
	pV = new double[dim*dim];
	memcpy(pV, eig.pV, dim*dim*sizeof(double));
	
	return *this;
}

ostream & operator << (ostream & os, const CIEig &eig)
{
	for ( INT i=0 ; i<eig.getDim() ; i++ )
		os << setw(16) << eig[i] << endl;
	return os;
}
