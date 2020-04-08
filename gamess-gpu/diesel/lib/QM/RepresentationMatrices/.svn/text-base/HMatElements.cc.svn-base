//***********************************************************************
//
//	Name:			HMatElements.cc
//
//	Description:	calculates Hamilton Matrix elements
//					
//
//	Author:			Michael Hanrath
//
//	Version:		0.0
//
//	Date:			10.09.1998
//
//
//
//
//
//***********************************************************************

#include "HMatElements.h"

#include "../Configuration/TableCase.h"

#include "../Math/SpinEigenFunctionDegeneration.h"
#include "RepDiag.h"

#include "../../Math/etc/Histogram.h"


#include <stdlib.h>

template <class MatrixType, class VectorType>
FourIndexIntegralContainer *
	HMatElements<MatrixType, VectorType>::int4Container = NULL;

template <class MatrixType, class VectorType>
TwoIndexIntegralContainer *
	HMatElements<MatrixType, VectorType>::int2Container = NULL;

template <class MatrixType, class VectorType>
MOType	HMatElements<MatrixType, VectorType>::maxMO = 0;


template <class MatrixType, class VectorType>
const SymmetryContainer *HMatElements<MatrixType, VectorType>::ijkk = NULL;
template <class MatrixType, class VectorType>
const SymmetryContainer *HMatElements<MatrixType, VectorType>::ijjk = NULL;
template <class MatrixType, class VectorType>
const SymmetryContainer *HMatElements<MatrixType, VectorType>::iijj = NULL;

template <class MatrixType, class VectorType>
double	HMatElements<MatrixType, VectorType>::core = 0;




template <class MatrixType, class VectorType>
void (HMatElements<MatrixType, VectorType>::*
	HMatElements<MatrixType, VectorType>::CbExMultTab[6][6])(
		CIVectors<VectorType> &, CIVectors<VectorType> &,
		INT, INT,
		IntegralType, IntegralType) const;


template <class MatrixType, class VectorType>
void (HMatElements<MatrixType, VectorType>::*
 HMatElements<MatrixType, VectorType>::CbMultTab[5])(
		CIVectors<VectorType> &, CIVectors<VectorType> &,
		INT, INT,
		IntegralType) const;



template <class MatrixType, class VectorType>
HMatElements<MatrixType, VectorType>::HMatElements(MOType _maxMO)
{
	maxMO = _maxMO;
	int4Container = NULL;
	int2Container = NULL;
	core = 0;
	this->CbMat = NULL;
	this->ExP3Mat = NULL;

}


template <class MatrixType, class VectorType>
HMatElements<MatrixType, VectorType>::HMatElements(
	FourIndexIntegralContainer * _int4Container,
	TwoIndexIntegralContainer * _int2Container,
	double _core)
{
//	printf("+++A %x\n", this);
//	fflush(stdout);
	int4Container = _int4Container;
	int2Container = _int2Container;
	core = _core;
	ijkk = (*int4Container)[FourIndexIntegralContainer::ijkk];
	ijjk = (*int4Container)[FourIndexIntegralContainer::ijjk];
	iijj = (*int4Container)[FourIndexIntegralContainer::iijj];
	this->CbMat = NULL;
	this->ExP3Mat = NULL;
/*	for ( INT i=0 ; i<5 ; i++ )
		for ( INT j=0 ; j<5 ; j++ )
			CbExMultTab[i][j] = (void (HMatElements<MatrixType, VectorType>::*)(
			CIVectors<CIVectorType> &, CIVectors<CIVectorType> &, INT, INT,
			IntegralType, IntegralType)) NULL;
*/			
	CbExMultTab[0][0] = &HMatElements<MatrixType, VectorType>::CbExMult_00;
	CbExMultTab[0][1] = &HMatElements<MatrixType, VectorType>::CbExMult_01;
	CbExMultTab[0][2] = &HMatElements<MatrixType, VectorType>::CbExMult_02;
	CbExMultTab[0][3] = &HMatElements<MatrixType, VectorType>::CbExMult_03;
	CbExMultTab[0][4] = &HMatElements<MatrixType, VectorType>::CbExMult_04;
	CbExMultTab[1][0] = &HMatElements<MatrixType, VectorType>::CbExMult_10;
	CbExMultTab[1][1] = &HMatElements<MatrixType, VectorType>::CbExMult_11;
	CbExMultTab[1][2] = &HMatElements<MatrixType, VectorType>::CbExMult_12;
	CbExMultTab[1][3] = &HMatElements<MatrixType, VectorType>::CbExMult_13;
	CbExMultTab[1][4] = &HMatElements<MatrixType, VectorType>::CbExMult_14;
	CbExMultTab[2][0] = &HMatElements<MatrixType, VectorType>::CbExMult_20;
	CbExMultTab[2][1] = &HMatElements<MatrixType, VectorType>::CbExMult_21;
	CbExMultTab[2][2] = &HMatElements<MatrixType, VectorType>::CbExMult_22;
	CbExMultTab[2][3] = &HMatElements<MatrixType, VectorType>::CbExMult_23;
	CbExMultTab[2][4] = &HMatElements<MatrixType, VectorType>::CbExMult_24;
	CbExMultTab[3][0] = &HMatElements<MatrixType, VectorType>::CbExMult_30;
	CbExMultTab[3][1] = &HMatElements<MatrixType, VectorType>::CbExMult_31;
	CbExMultTab[3][2] = &HMatElements<MatrixType, VectorType>::CbExMult_32;
	CbExMultTab[3][3] = &HMatElements<MatrixType, VectorType>::CbExMult_33;
	CbExMultTab[3][4] = &HMatElements<MatrixType, VectorType>::CbExMult_34;
	CbExMultTab[4][0] = &HMatElements<MatrixType, VectorType>::CbExMult_40;
	CbExMultTab[4][1] = &HMatElements<MatrixType, VectorType>::CbExMult_41;
	CbExMultTab[4][2] = &HMatElements<MatrixType, VectorType>::CbExMult_42;
	CbExMultTab[4][3] = &HMatElements<MatrixType, VectorType>::CbExMult_43;
	CbExMultTab[4][4] = &HMatElements<MatrixType, VectorType>::CbExMult_44;
	CbExMultTab[5][5] = &HMatElements<MatrixType, VectorType>::CbExMult_55;

	CbMultTab[0] = &HMatElements<MatrixType, VectorType>::CbMult_0;
	CbMultTab[1] = &HMatElements<MatrixType, VectorType>::CbMult_1;
	CbMultTab[2] = &HMatElements<MatrixType, VectorType>::CbMult_2;
	CbMultTab[3] = &HMatElements<MatrixType, VectorType>::CbMult_3;
	CbMultTab[4] = &HMatElements<MatrixType, VectorType>::CbMult_4;
}

template <class MatrixType, class VectorType>
HMatElements<MatrixType, VectorType>::~HMatElements()
{
}


#include "HMatElements.Dens.cch"
#include "HMatElements.Diag.cch"
#include "HMatElements.Mult.cch"
#include "HMatElements.Get.cch"
#include "HMatElements.PT.cch"


template class HMatElements<float, float>;
template class HMatElements<double, float>;
template class HMatElements<float, double>;
template class HMatElements<double, double>;
