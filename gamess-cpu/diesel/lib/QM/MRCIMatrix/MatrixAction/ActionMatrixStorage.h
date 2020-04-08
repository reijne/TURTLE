//***********************************************************************
//
//	Name:			ActionMatrixStorage.h
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


#ifndef __ActionMatrixStorage_h
#define __ActionMatrixStorage_h

#include "../../../../config.h"


#include "../../MRTree/MatrixAction.h"
#include "../../RepresentationMatrices/HMatElements.h"
#include "../MatrixStorage.h"
#include "../../Configuration/DiffConf.h"
#include "../../Configuration/Configuration.h"

template <class MatrixType, class VectorType>
class ActionMatrixStorage : 
	public MatrixAction<HMatElements<MatrixType, VectorType> > {

public:
	ActionMatrixStorage(MatrixStorage<MatrixType, VectorType> *stored) :
		stored(stored) {}
		
	void	doit(
		const HMatElements<MatrixType, VectorType> & repMats,
		INT safStartA, INT safStartB,
		const DiffConf<MOType> & diffConf, const TableCase<MOType> & tablecase)
	{
	MatrixType*  p = new MatrixType[repMats.getNumberOfRows()*repMats.getNumberOfColumns()];

		repMats.getMatrix(p, diffConf, tablecase);


	MatrixType	*pp = p;
		for ( INT k=0 ; k<repMats.getNumberOfRows() ; k++ )
			for ( INT l=0 ; l<repMats.getNumberOfColumns() ; l++ )
			{
				if ( tablecase.getP()!=5 || k>=l )
					stored->addEntry(safStartA+k, safStartB+l, *pp);
				pp++;
			}
	delete[] p;
	}
MatrixStorage<MatrixType, VectorType> *stored;
};





#endif
