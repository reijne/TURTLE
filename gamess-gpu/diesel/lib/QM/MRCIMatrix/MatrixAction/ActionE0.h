//***********************************************************************
//
//	Name:			ActionE0.h
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


#ifndef __ActionE0_h
#define __ActionE0_h

#include "../../../../config.h"


#include "../../MRTree/MatrixAction.h"
#include "../../RepresentationMatrices/MRMPH0MatElements.h"
#include "../../Configuration/DiffConf.h"
#include "../../Configuration/Configuration.h"


template <class MatrixType, class VectorType>
class ActionE0 : 
	public MatrixAction<MRMPH0MatElements<MatrixType, VectorType> > {

public:
	ActionE0(const VectorType	*Psi0EV, MatrixType	&E0) :
		Psi0EV(Psi0EV), E0(E0) {}
		
	void	doit(
		const MRMPH0MatElements<MatrixType, VectorType> & repMats,
		INT safStartA, INT safStartB,
		const DiffConf<MOType> & diffConf, const TableCase<MOType> & tablecase)
	{
	MatrixType* p = new MatrixType[repMats.getNumberOfRows()*repMats.getNumberOfColumns()];

		repMats.getMatrix(p, diffConf, tablecase);

	const MatrixType	*pp = p;
		for ( INT i=0 ; i<repMats.getNumberOfRows() ; i++ )
			for ( INT j=0 ; j<repMats.getNumberOfColumns() ; j++ )
				E0 += *pp++ * Psi0EV[safStartA+i] * Psi0EV[safStartB+j];
	delete[] p;
	}
const VectorType	*Psi0EV;
MatrixType	&E0;
};



#endif
