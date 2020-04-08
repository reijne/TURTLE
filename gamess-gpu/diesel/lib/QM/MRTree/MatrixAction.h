//***********************************************************************
//
//	Name:			MatrixAction.h
//
//	Description:	
//
//	Author:			Michael Hanrath
//
//	Version:		0.0
//
//	Date:			06.10.1998
//
//
//
//
//
//***********************************************************************


#ifndef __MatrixAction_h
#define __MatrixAction_h

#include "../../../config.h"


template <class MatrixElementCalculator> 
class MatrixAction {
public:
	MatrixAction() {}
	virtual ~MatrixAction() {}
	
	
	virtual void	doit(
		const MatrixElementCalculator & repMats,
		INT safStartA, INT safStartB,
		const DiffConf<MOType> &, const TableCase<MOType> &) = 0;

};





#endif
