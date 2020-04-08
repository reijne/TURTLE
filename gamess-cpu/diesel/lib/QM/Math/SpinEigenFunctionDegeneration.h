//***********************************************************************
//
//	Name:	SpinEigenFunctionDegeneration.h
//
//	Description:	
//
//	Author:	Michael Hanrath
//
//	Version:	0.0
//
//	Date:	23.08.1996
//
//
//
//
//***********************************************************************

#ifndef __SPINEIGENFUNCTIONDEGENERATION_H
#define __SPINEIGENFUNCTIONDEGENERATION_H

#include "../../../config.h"

#include <iostream>
using std::ostream;

class SpinEigenFunctionDegeneration {
public:
	SpinEigenFunctionDegeneration(INT multiplicity, INT maxOpenShells);

	SpinEigenFunctionDegeneration(const SpinEigenFunctionDegeneration &);

	~SpinEigenFunctionDegeneration();


	INT	operator () (INT openShells) const;
	
	double	getTotalSpin() const;
	INT	getMultiplicity() const;
	INT	getMaxOpenShells() const;

	const INT *getP() const;

//----------------------------------------------------------------------	
	
	friend ostream& operator<<(ostream& s,
		const SpinEigenFunctionDegeneration & seigs);
	
//----------------------------------------------------------------------	

	
private:
	INT	fak(INT);
	INT	nuebk(INT n, INT k);
	
INT	totalSpin2;
INT	maxOpenShells;
INT	*p;
};



inline
INT	SpinEigenFunctionDegeneration::operator() (INT openShells) const
{	return	p[openShells >> 1];	}


inline
double	SpinEigenFunctionDegeneration::getTotalSpin() const
{	return	totalSpin2 / 2.0;	}

inline
INT	SpinEigenFunctionDegeneration::getMultiplicity() const
{	return	totalSpin2 + 1;	}

inline
INT	SpinEigenFunctionDegeneration::getMaxOpenShells() const
{	return	maxOpenShells;	}


inline
const INT *SpinEigenFunctionDegeneration::getP() const
{	return	p;	}



#endif
