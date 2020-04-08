//***********************************************************************
//
//	Name:			Roots.cc
//
//	Description:	
//
//	Author:			Michael Hanrath
//
//	Version:		0.0
//
//	Date:			11.02.1997
//
//***********************************************************************


#include "Roots.h"


#include <string>

#include <iomanip>

using namespace std;

Roots::Roots(INT _n)
{
	n = _n;
	roots = new TRoot[n];
	memset(roots, 0, n*sizeof(TRoot));
	for ( INT i=0 ; i<n ; i++ )
	{
		roots[i].rootNumber = i;
		roots[i].status = Uninitialized;
		roots[i].root = 0;
		roots[i].refRoot = 0;
	}
}

Roots::Roots(INT _n, INT rootNumber1 ...)
{
	n = _n;
	roots = new TRoot[n];
	memset(roots, 0, n*sizeof(TRoot));
	
va_list	ap;
	va_start(ap, rootNumber1);
	roots[0].rootNumber = rootNumber1;
	for ( INT i=0 ; i<n ; i++ )
	{
		if ( i )
			roots[i].rootNumber = va_arg(ap, INT);
		roots[i].status = Uninitialized;
		roots[i].root = 0;
		roots[i].refRoot = 0;
	}
}

Roots::Roots(INT _n, const INT *rootNumbers)
{
	n = _n;
	roots = new TRoot[n];
	memset(roots, 0, n*sizeof(TRoot));
	for ( INT i=0 ; i<n ; i++ )
	{
		roots[i].rootNumber = rootNumbers[i] - 1;
		roots[i].status = Uninitialized;
		roots[i].root = 0;
		roots[i].refRoot = 0;
	}
}



Roots::~Roots()
{
	delete	roots;
}


char	*Roots::getConvergenceName(ConvergenceStatus status) const
{
static char *names[3] = {	"uninitialized", "converged", "not converged" };
	return names[status];
}


ostream& operator<<(ostream & s, const Roots &r)
{
	s << " no.  root   convergence status             value" << endl;
	s << "---------------------------------------------------" << endl;
	for ( INT i=0 ; i<r.getNumberOfRoots() ; i++ )
	{
		s << setw(4) << i+1;
		s << setw(6) << r.getRootNumber(i)+1;
		s << setw(21) << r.getConvergenceName(r.getConvergenceStatus(i));
		if ( r.getConvergenceStatus(i)!=Roots::Uninitialized )
			s << setw(18) << r.getRoot(i);
		else
			s << setw(18) << "---";
		s << endl;
	}
	return s;
}

