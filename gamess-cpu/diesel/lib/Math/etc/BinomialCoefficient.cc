//***********************************************************************
//
//	Name:	BinomialCoefficient.cc
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


#include "BinomialCoefficient.h"

#include <string>

using namespace std; 

BinomialCoefficient::BinomialCoefficient(INT _maxN)
{
	maxN = _maxN;
	p = new INT [maxN*maxN];

	memset(p, 0, maxN*maxN*sizeof(INT));
	for ( INT i=0 ; i<maxN ; i++ )
	{	set(i, 0, 1);
		set(i, 1, i);
	}
	for ( INT i=2 ; i<maxN ; i++ )
		for ( INT j=2 ; j<=i ; j++ )
			set(i, j, get(i-1, j) + get(i-1, j-1));
}



BinomialCoefficient::~BinomialCoefficient()
{
	delete[] p;
}
