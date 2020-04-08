//***********************************************************************
//
//	Name:	BinomialCoefficient.h
//
//	Description:	
//
//	Author:	Michael Hanrath
//
//	Version:	0.0
//
//	Date:	07.08.1996
//
//
//
//
//***********************************************************************

#ifndef __BINOMIALCOEFFICIENT_H
#define __BINOMIALCOEFFICIENT_H

#include "../../../config.h"



class BinomialCoefficient {
public:
	BinomialCoefficient(INT maxN);
	~BinomialCoefficient();
	
//----------------------------------------------------------------------	
	
	INT	get(INT n, INT k) const;
	
//----------------------------------------------------------------------	
	
	
private:
	void	set(INT n, INT k, INT noverk);

INT	*p;
INT	maxN;
};



inline
INT	BinomialCoefficient::get(INT n, INT k) const
{	return p[n*maxN + k];	}

inline
void	BinomialCoefficient::set(INT n, INT k, INT noverk)
{	p[n*maxN + k] = noverk;	}



#endif
