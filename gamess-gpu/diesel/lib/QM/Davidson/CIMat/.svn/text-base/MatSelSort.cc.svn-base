//***********************************************************************
//
//	Name:			MatSelSort.cc
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

#include "MatSelSort.h"

#include <stdlib.h>
#include <string>
#include <iostream>
#include <iomanip>

using namespace std;

MatSelSort::MatSelSort(INT dim) :
	MatSel(dim)
{
	e = new double[dim];
	memset(e, 0, dim*sizeof(double));
}


MatSelSort::~MatSelSort()
{
	if ( e )
		delete e;
}

int	static cmp(const void *_p1, const void *_p2)
{
const MatSelSort::TRecord *p1 = (const MatSelSort::TRecord *) _p1;
const MatSelSort::TRecord *p2 = (const MatSelSort::TRecord *) _p2;
	if ( p1->v>p2->v )
		return 1;
	if ( p1->v<p2->v )
		return -1;
	return 0;
}

void	MatSelSort::sort()
{
TRecord	*q = new TRecord[dim];

	for ( INT i=0 ; i<dim ; i++ )
	{
		q[i].n = p[i];
		q[i].v = e[i];
	}

	qsort(q, dim, sizeof(TRecord), cmp);

	for ( INT i=0 ; i<dim ; i++ )
	{
		p[i] = q[i].n;
		e[i] = q[i].v;
	}

	delete q;	
}


ostream & operator << (ostream & os, const MatSelSort & sel)
{
	for ( INT i=0 ; i<sel.getDim() ; i++ )
		os << i << " --> " << sel[i] << " (" << sel.getValue(i) << ")" << endl;
	return os;
}
