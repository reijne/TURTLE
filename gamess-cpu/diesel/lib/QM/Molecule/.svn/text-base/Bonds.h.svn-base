//***********************************************************************
//
//	Name:			Bonds.h
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

#ifndef __BONDS_H
#define __BONDS_H

#include "../../../config.h" 
#include "../../Container/SLList.h"

#include <iostream>
using std::ostream;

class Molecule;
class Atom;


// FD class ostream; 

struct	TBond {
	TBond(double _dist, INT _order, const Atom *_a1, const Atom *_a2) : 
		dist(_dist), order(_order), a1(_a1), a2(_a2) {}

	double	dist;
	INT	order;
	const Atom	*a1, *a2;
	};

class Bonds : public SLList<TBond *> {
public:

	Bonds(const Molecule &);
	~Bonds();
	

	friend ostream & operator << (ostream &, const Bonds &);

private:

};






#endif
