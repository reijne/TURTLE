//***********************************************************************
//
//	Name:			ActionMultiplication.h
//
//	Description:	
//
//	Author:			Michael Hanrath
//
//	Version:		0.0
//
//	Date:			07.10.1998
//
//
//
//
//
//***********************************************************************


#ifndef __ActionMultiplication_h
#define __ActionMultiplication_h

#include "../../../../config.h"

#include "../../RepresentationMatrices/CIVectors.h"

#include "../../MRTree/MatrixAction.h"
#include "../../RepresentationMatrices/HMatElements.h"
#include "../../Configuration/DiffConf.h"
#include "../../Configuration/Configuration.h"

template <class MatrixType, class VectorType>
class ActionMultiplication : 
	public MatrixAction<HMatElements<MatrixType, VectorType> > {

public:
	ActionMultiplication(CIVectors<VectorType> &CIx, CIVectors<VectorType> &CIy) :
		CIx(CIx), CIy(CIy) {}
		
	void	doit(
		const HMatElements<MatrixType, VectorType> & repMats,
		INT safStartA, INT safStartB,
		const DiffConf<MOType> & diffConf, const TableCase<MOType> & tablecase)
	{
	repMats.Mult(CIy, CIx, 
		safStartA, safStartB, 
		diffConf, tablecase);
	}
CIVectors<VectorType> &CIx;
CIVectors<VectorType> &CIy;
};




#endif
