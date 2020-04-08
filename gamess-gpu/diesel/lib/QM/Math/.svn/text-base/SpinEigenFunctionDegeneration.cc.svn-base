//***********************************************************************
//
//	Name:	SpinEigenFunctionDegeneration.cc
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

#include "SpinEigenFunctionDegeneration.h"

#include <math.h>
#include <string>


SpinEigenFunctionDegeneration::SpinEigenFunctionDegeneration
	(INT _multiplicity, INT _maxOpenShells)
{
	totalSpin2 = _multiplicity - 1;
	maxOpenShells = _maxOpenShells;
	p = new INT [maxOpenShells/2+1];
	memset(p, 0, (maxOpenShells/2+1)*sizeof(INT));
	for ( INT i=totalSpin2 ; i<=maxOpenShells ; i+=2 )
		p[i/2] = nuebk(i, i/2 - totalSpin2/2) - nuebk(i, i/2 - totalSpin2/2 - 1);
		
/*		(totalSpin2+1) * fak(i) /
			(fak(i/2 + totalSpin2/2 + 1 + (i&1) ) * fak(i/2 - totalSpin2/2));
*/
}


SpinEigenFunctionDegeneration::SpinEigenFunctionDegeneration(
	const SpinEigenFunctionDegeneration &s)
{
	totalSpin2 = s.totalSpin2;
	maxOpenShells = s.maxOpenShells;
	p = new INT [maxOpenShells/2+1];
	memcpy(p, s.p, (maxOpenShells/2+1)*sizeof(INT));
}



SpinEigenFunctionDegeneration::~SpinEigenFunctionDegeneration()
{
	delete[] p;
}

INT	SpinEigenFunctionDegeneration::fak(INT i)
{	return (i>1) ? i*fak(i-1) : 1;	}


INT	SpinEigenFunctionDegeneration::nuebk(INT n, INT k)
{
	if ( k<0 )
		return 0;
INT	p = 1;
INT	j = 1;
	for ( INT i=n-k+1 ; i<=n ; i++ )
	{	p *= i;
		if ( j<=k )
			p /= j++;
	}
	for ( ; j<=k ; j++ )
		p /= j;
	return p;
}

ostream& operator<<(ostream& s,	const SpinEigenFunctionDegeneration & seigs)
{
	s << "multiplicity: " << seigs.getMultiplicity() << std::endl;
	s << "total spin  : " << seigs.getTotalSpin() << std::endl;
	s << "open shells                      : ";
	for ( INT i=seigs.getMultiplicity()-1 ; i<=seigs.getMaxOpenShells() ; i+=2 )
		s << i << " ";
	s << std::endl;
	s << "degenaracy of spin eigenfunctions: "; 
	for ( INT i=seigs.getMultiplicity()-1 ; i<=seigs.getMaxOpenShells() ; i+=2 )
		s << seigs(i) << " ";
	s << std::endl;
	return s;
}
