//***********************************************************************
//
//	Name:			Roots.h
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

#ifndef __ROOTS_H
#define __ROOTS_H

#include <iostream>
using std::istream;
using std::ostream;

#include <stdarg.h>

#include "../../../config.h"

class Roots {
public:
	Roots(INT n);
	Roots(INT n, INT rootNumber1 ...);
	Roots(INT n, const INT *rootNumbers);
	
	~Roots();

enum	ConvergenceStatus {
	Uninitialized, Converged, NotConverged };

//--------------------------------------------------------------------------

	INT	getNumberOfRoots() const;
	INT	getRootNumber(INT i) const;
	ConvergenceStatus	getConvergenceStatus(INT i) const;
	double	getRoot(INT i) const;
	double	getRefRoot(INT i) const;

	char	*getConvergenceName(ConvergenceStatus status) const;

	void	setRootNumber(INT i, INT n) const;
	void	setConvergenceStatus(INT i, ConvergenceStatus status);
	void	setRoot(INT i, double root);
	void	setRefRoot(INT i, double root);

//--------------------------------------------------------------------------

	friend ostream& operator<<(ostream & s, const Roots &);
	friend istream& operator>>(istream & s, Roots &);

//--------------------------------------------------------------------------

private:
INT	n;							//	number of roots
struct TRoot {
	INT	rootNumber;				//	number of root
	ConvergenceStatus	status;	//	status of convergence
	double	root;				//	eigenvalue (root) (in iteration process)
	double	refRoot;			//	eigenvalue (in reference space)
	} *roots;
	
	

};


inline
INT	Roots::getNumberOfRoots() const
{	return	n;	}

inline
INT	Roots::getRootNumber(INT i) const
{	return	roots[i].rootNumber;	}

inline
void	Roots::setRootNumber(INT i, INT n) const
{	roots[i].rootNumber = n;	}

inline
Roots::ConvergenceStatus	Roots::getConvergenceStatus(INT i) const
{	return	roots[i].status;	}

inline
double	Roots::getRoot(INT i) const
{	return	roots[i].root;	}

inline
double	Roots::getRefRoot(INT i) const
{	return	roots[i].refRoot;	}


inline
void	Roots::setConvergenceStatus(INT i, ConvergenceStatus status)
{	roots[i].status = status;	}

inline
void	Roots::setRoot(INT i, double root)
{	roots[i].root = root;	}

inline
void	Roots::setRefRoot(INT i, double root)
{	roots[i].refRoot = root;	}

#endif
