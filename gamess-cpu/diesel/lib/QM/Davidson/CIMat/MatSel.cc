//***********************************************************************
//
//	Name:			MatSel.cc
//
//	Description:	
//
//	Author:			Michael Hanrath
//
//	Version:		0.0
//
//	Date:			01. Aug 1998
//
//***********************************************************************

#include "MatSel.h"

#include <string>
#include <iostream>
#include <iomanip>

using namespace std;

MatSel::MatSel(INT _dim)
{
	dim = _dim;
	p = new INT[dim];
	neutral();
}

MatSel::~MatSel()
{
	if ( p )
		delete p;
}

INT	MatSel::del(INT where)
{
INT	*pp = new INT[dim-1];
INT	v = 0;
INT	j = 0;

	for ( INT i=0 ; i<dim ; i++ )
		if ( i!=where )
			pp[j++] = p[i];
		else
			v = p[i];
			
	delete p;
	dim--;
	p = pp;
	return v;
}


void	MatSel::ins(INT where, INT what)
{
INT	*pp = new INT[dim+1];
INT	j = 0;

	for ( INT i=0 ; i<dim+1 ; i++ )
		if ( i!=where )
			pp[i] = p[j++];
		else
			pp[i] = what;
			
	delete p;
	dim++;
	p = pp;
}



void	MatSel::neutral()
{
	for ( INT i=0 ; i<dim ; i++ )
		p[i] = i;
}



ostream & operator << (ostream & os, const MatSel & sel)
{
	for ( INT i=0 ; i<sel.dim ; i++ )
		os << i << " --> " << ((MatSel &) sel)[i] << endl;
	return os;
}
