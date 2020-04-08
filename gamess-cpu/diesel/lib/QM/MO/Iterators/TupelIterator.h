//***********************************************************************
//
//	Name:			TupelIterator.h
//
//	Description:	iterators for MO lists
//
//	Author:			Michael Hanrath
//
//	Version:		0.0
//
//	Date:			23.07.1997
//
//
//
//
//
//***********************************************************************

#ifndef __TupelIterator_h
#define __TupelIterator_h

#include "../../../../config.h"

#include <iostream>
using std::ostream;

class TupelIterator {
public:
	TupelIterator(INT nElectrons);
	~TupelIterator() {}


	void	first();
	void	next();
	INT	isEnd() const;
	
	
	INT	getNumberOfOpenShells() const;
	INT	getNumberOfClosedShells() const;
	INT	getNumberOfShells() const;

	//----------------------------------------------------------------	

	friend ostream& operator<<(ostream & s, const TupelIterator &);

	//----------------------------------------------------------------

private:
INT	n;			// number of electrons
INT	iter;		// iteration variable
};

inline
TupelIterator::TupelIterator(INT nElectrons)
{
	n = nElectrons;
}

inline
void	TupelIterator::first()
{	iter = n;	}

inline
void	TupelIterator::next()
{	iter -= 2;	}

inline
INT	TupelIterator::isEnd() const
{	return iter<0;	}



inline
INT	TupelIterator::getNumberOfOpenShells() const
{	return	iter;	}

inline
INT	TupelIterator::getNumberOfClosedShells() const
{	return	(n-iter)/2;	}

inline
INT	TupelIterator::getNumberOfShells() const
{	return	(n+iter)/2;	}



#endif


