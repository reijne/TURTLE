//***********************************************************************
//
//	Name:			MRMPH0MatElements.cc
//
//	Description:	calculates Multi Reference Moller-Plesset
//					Hamilton Matrix elements
//					
//
//	Author:			Michael Hanrath
//
//	Version:		0.0
//
//	Date:			14.09.1998
//
//
//
//***********************************************************************




#include "MRMPH0MatElements.h"
#include "HMatElements.h"

#include "../IntegralContainer/MRFockMatrix.h"




template <class MatrixType, class VectorType>
TwoIndexIntegralContainer *
	MRMPH0MatElements<MatrixType, VectorType>::int2Container = NULL;

template <class MatrixType, class VectorType>
const SymmetryContainer *
	MRMPH0MatElements<MatrixType, VectorType>::iijj = NULL;

template <class MatrixType, class VectorType>
FourIndexIntegralContainer *
	MRMPH0MatElements<MatrixType, VectorType>::int4Container = NULL;

template <class MatrixType, class VectorType>
double
	MRMPH0MatElements<MatrixType, VectorType>::core = 0;

template <class MatrixType, class VectorType>
const MRFockMatrix<VectorType>	*
	MRMPH0MatElements<MatrixType, VectorType>::mrFockMatrix = NULL;

template <class MatrixType, class VectorType>
MatrixType
	MRMPH0MatElements<MatrixType, VectorType>::E0 = 0;


template <class MatrixType, class VectorType>
MRMPH0MatElements<MatrixType, VectorType>::MRMPH0MatElements(
	const MRFockMatrix<VectorType> * _mrFockMatrix,
	TwoIndexIntegralContainer * _int2Container,
	FourIndexIntegralContainer * _int4Container,
	double _core,
	const MatrixType _E0)
{
	this->occupiedMemory = 0;
	mrFockMatrix = _mrFockMatrix;
	int2Container = _int2Container;
	int4Container = _int4Container;
	core = _core;
	E0 = _E0;
	iijj = (*int4Container)[FourIndexIntegralContainer::iijj];
}





template <class MatrixType, class VectorType>
MRMPH0MatElements<MatrixType, VectorType>::~MRMPH0MatElements()
{
}




template <class MatrixType, class VectorType>
MatrixType	MRMPH0MatElements<MatrixType, VectorType>::getDiag(
	const Configuration<MOType> & same)
{
MatrixType	diag = core;
//	cout << "diag1 " << diag <<endl;
	for ( INT i=0 ; i<same.getNumberOfOpenShells() ; i++ )
	{
MOType	moi;
		moi = same.getOpenShell(i);
		diag += (*mrFockMatrix)(moi, moi);
	}
//	cout << "diag2 " << diag <<endl;
	for ( INT k=0 ; k<same.getNumberOfClosedShells() ; k++ )
	{
MOType	mok;
		mok = same.getClosedShell(k);
		diag += 2*(*mrFockMatrix)(mok, mok);
	}
//	cout << "diag3 " << diag <<endl;

	return diag;
}




#include "MRMPH0MatElements.Get.cch"
#include "MRMPH0MatElements.Mult.cch"




template class MRMPH0MatElements<float, float>;
template class MRMPH0MatElements<double, float>;
template class MRMPH0MatElements<float, double>;
template class MRMPH0MatElements<double, double>;
