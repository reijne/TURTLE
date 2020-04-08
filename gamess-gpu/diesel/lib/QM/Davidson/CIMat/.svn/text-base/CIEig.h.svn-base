#include "../../../../config.h"
//***********************************************************************
//
//	Name:			CIEig.h
//
//	Description:	
//
//	Author:			Michael Hanrath
//
//	Version:		0.0
//
//	Date:			15. Aug 1998
//
//***********************************************************************

#ifndef __CIEig_h
#define __CIEig_h

#include <iostream>
using std::ostream;

/* FD class ostream; */
class CIMat;


class CIEig {
public:
enum Method { FullMethod, DavidsonMethod };
	CIEig(const CIMat &);
	CIEig(const CIMat &, INT nRoots);
	~CIEig();

	CIEig(const CIEig &);
	CIEig & operator = (const CIEig &);

	INT	getDim() const;


	double	operator [] (INT) const;
	
	double	getEV(INT n, INT component) const;


    	friend ostream & operator << (ostream &, const CIEig &);

private:

INT	nRoots;
INT	dim;
double	*pE;
double	*pV;

};



inline
INT	CIEig::getDim() const
{	return	dim;	}


inline
double	CIEig::operator [] (INT i) const
{	return	pE[i];	}


inline
double	CIEig::getEV(INT n, INT component) const
{	return	pV[n*dim + component];	}


#endif
