//***********************************************************************
//
//	Name:			ActionDensity.h
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


#ifndef __ActionDensity_h
#define __ActionDensity_h

#include "../../../../config.h"


#include "../../MRTree/MatrixAction.h"
#include "../../RepresentationMatrices/HMatElements.h"
#include "../../Configuration/DiffConf.h"
#include "../../Configuration/Configuration.h"


template <class MatrixType, class VectorType>
class ActionDensity : 
	public MatrixAction<HMatElements<MatrixType, VectorType> > {

public:
	ActionDensity(CIVectors<VectorType> &CIx, CIVectors<VectorType> &CIy,
		VectorType	**densMat) :
		CIx(CIx), CIy(CIy), densMat(densMat) {}
		
	void	doit(
		const HMatElements<MatrixType, VectorType> & repMats,
		INT safStartA, INT safStartB,
		const DiffConf<MOType> & diffConf, const TableCase<MOType> & tablecase)
	{
		repMats.Dens(&CIy, &CIx, 
			safStartA, safStartB, 
			diffConf, tablecase, densMat);
	}
CIVectors<VectorType> &CIx;
CIVectors<VectorType> &CIy;
VectorType	**densMat;
};





#endif
