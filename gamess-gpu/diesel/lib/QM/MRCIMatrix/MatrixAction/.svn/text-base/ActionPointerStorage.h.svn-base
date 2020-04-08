//***********************************************************************
//
//	Name:			ActionPointerStorage.h
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


#ifndef __ActionPointerStorage_h
#define __ActionPointerStorage_h

#include "../../../../config.h"


#include "../../MRTree/MatrixAction.h"
#include "../../RepresentationMatrices/HMatElements.h"
#include "../../Configuration/DiffConf.h"
#include "../../Configuration/Configuration.h"



template <class MatrixType, class VectorType>
class ActionPointerStorage : 
	public MatrixAction<HMatElements<MatrixType, VectorType> > {

public:
	ActionPointerStorage(MatrixType	*HMat, INT SAFs) :
		HMat(HMat), SAFs(SAFs) {}
		
	void	doit(
		const HMatElements<MatrixType, VectorType> & repMats,
		INT safStartA, INT safStartB,
		const DiffConf<MOType> & diffConf, const TableCase<MOType> & tablecase)
	{
	MatrixType* p = new MatrixType[repMats.getNumberOfRows()*repMats.getNumberOfColumns()];

		repMats.getMatrix(p, diffConf, tablecase);

	MatrixType	*pp = p;
		for ( INT k=0 ; k<repMats.getNumberOfRows() ; k++ )
			for ( INT l=0 ; l<repMats.getNumberOfColumns() ; l++ )
				HMat[(safStartA+k)*SAFs + safStartB+l] = 
					HMat[(safStartB+l)*SAFs + safStartA+k] = *pp++;
	delete[] p;
	}
MatrixType	*HMat;								// pointer to hamilton matrix
INT	SAFs;
};



#endif
