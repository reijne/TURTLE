//***********************************************************************
//
//	Name:			Bonds.cc
//
//	Description:	
//
//	Author:			Michael Hanrath
//
//	Version:		0.0
//
//	Date:			12. Feb 1998
//
//***********************************************************************

#include "Bonds.h"
#include "Molecule.h"
#include "Atom.h"

#include <iostream>

using namespace std;

Bonds::Bonds(const Molecule & m)
{
	for ( INT i=0 ; i<m.getNumberOfAtoms() ; i++ )
	{
		for ( INT j=i+1 ; j<m.getNumberOfAtoms() ; j++ )
		{
		TBond	bond = m[i].getBond(m[j]);
			if ( bond.order )
			{
			TBond	*Pbond = new TBond(bond);
				append(Pbond);
			}
		}
	}
		
}



Bonds::~Bonds()
{
Pix	pix = first();
	while ( pix )
	{
		delete (*this)(pix);
		next(pix);
	}
}


ostream & operator << (ostream & s, const Bonds & b)
{
Pix	pix = b.first();
	while ( pix )
	{
	TBond *p = b(pix);
		s << *p->a1 << " <--> " << *p->a2 << ", dist=" << p->dist <<
			", order=" << p->order << endl;
		b.next(pix);
	}

	return s;
}


#include "../../Container/SLList.cc"
template class SLList<TBond *>;

