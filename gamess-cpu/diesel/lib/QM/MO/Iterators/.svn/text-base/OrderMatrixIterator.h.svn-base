//***********************************************************************
//
//	Name:			OrderMatrixIterator.h
//
//	Description:	iterator for order relation matrices
//
//	Author:			Michael Hanrath
//
//	Version:		0.0
//
//	Date:			22.07.1997
//
//
//
//
//
//***********************************************************************

#ifndef __OrderMatrixIterator_h
#define __OrderMatrixIterator_h

#include "../../../../config.h"

#include "../MOType.h"
#include "../GeneralizedMO.h"
#include "../../Configuration/Configuration.h"
#include "../../Configuration/DiffConf.h"
#include "MOListIterator.h"

#include <iostream>

class OrderMatrixIterator {
public:
	OrderMatrixIterator(INT rows, INT columns);
	~OrderMatrixIterator();

	//----------------------------------------------------------------	

	INT	getNumberOfRows() const;
	INT	getNumberOfColumns() const;
	INT	getDim() const;
	
	INT	operator() (INT row, INT column) const;

	INT	isEqualInColumn(INT i) const;
	INT	isEqualInRow(INT i) const;

	//----------------------------------------------------------------	

	INT	calcMOListIterator(
		MOType start, MOType end, 
		INT tupelClosed,
		const MOType *baseMOsInIrRep,
		MOListIterator & molist,
		Configuration<MOType> &extNotRunning,
		DiffConf<GeneralizedMO> &dcGen,
		INT &sig,
		IrRep irrep) const;
		
	//----------------------------------------------------------------	

	void	first();
	void	next();
	INT	isEnd() const;

	//----------------------------------------------------------------	

	friend ostream& operator<<(ostream & s, const OrderMatrixIterator &);

	//----------------------------------------------------------------
private:
INT	rows;		// number of rows
INT	cols;		// number of columns
INT	dim;		// dimension of matrix (r*c)
INT	*p;			// relation matrix
INT	*row0;		// number of zeros in specific row
INT	*col0;		// number of zeros in specific column

INT	i;			// loop variable
INT	end;		// flag if end
};


inline
INT	OrderMatrixIterator::isEnd() const
{	return end==2;	}


inline
INT	OrderMatrixIterator::getNumberOfRows() const
{	return	rows;	}

inline
INT	OrderMatrixIterator::getNumberOfColumns() const
{	return	cols;	}
	
inline
INT	OrderMatrixIterator::getDim() const
{	return	dim;	}

inline
INT	OrderMatrixIterator::operator() (INT row, INT column) const
{	return p[row*cols + column];	}

inline
INT	OrderMatrixIterator::isEqualInColumn(INT i) const
{	return	col0[i]>0;	}

inline
INT	OrderMatrixIterator::isEqualInRow(INT i) const
{	return	row0[i]>0;	}


#endif
