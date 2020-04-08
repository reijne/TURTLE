//***********************************************************************
//
//	Name:			ActionE1.h
//
//	Description:	
//
//	Author:			Michael Hanrath
//
//	Version:		0.0
//
//	Date:			12.10.1998
//
//
//
//
//
//***********************************************************************


#ifndef __ActionE1_h
#define __ActionE1_h

#include "../../../../config.h"


#include "../../MRTree/MatrixAction.h"
#include "../../RepresentationMatrices/HMatElements.h"
#include "../../Configuration/DiffConf.h"
#include "../../Configuration/Configuration.h"


template <class MatrixType, class VectorType>
class ActionE1 : 
	public MatrixAction<HMatElements<MatrixType, VectorType> > {

public:
	ActionE1(const VectorType	*Psi0EV, MatrixType	&E1) :
		Psi0EV(Psi0EV), E1(E1) {}
		
	void	doit(
		const HMatElements<MatrixType, VectorType> & repMats,
		INT safStartA, INT safStartB,
		const DiffConf<MOType> & diffConf, const TableCase<MOType> & tablecase)
	{
	MatrixType* p = new MatrixType[repMats.getNumberOfRows()*repMats.getNumberOfColumns()];

		repMats.getMatrix(p, diffConf, tablecase);

	const MatrixType	*pp = p;
		for ( INT i=0 ; i<repMats.getNumberOfRows() ; i++ )
			for ( INT j=0 ; j<repMats.getNumberOfColumns() ; j++ )
				E1 += *pp++ * Psi0EV[safStartA+i] * Psi0EV[safStartB+j];
	delete[] p;
	}
const VectorType	*Psi0EV;
MatrixType	&E1;
};



#endif
