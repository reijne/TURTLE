//***********************************************************************
//
//	Name:			Atom.h
//
//	Description:	
//
//	Author:			Michael Hanrath
//
//	Version:		0.0
//
//	Date:			07. Feb 1998
//
//***********************************************************************

#ifndef _ATOM_H
#define _ATOM_H

#include "../../../config.h"

#include "Point.h"
#include "TypeClass.h"

#include "Bonds.h"

#include <iostream>

#include "../../Container/String.h"

//FD class	ostream;

class Atom : public TypeClass, public Point  {
public:
	Atom();
	Atom(istream &s);
	Atom(String label, double charge, double x, double y, double z);


	TBond	getBond(const Atom &) const;
	INT getNr() const { return nr; }

	INT	operator == (const Atom & a) const;
	INT	operator != (const Atom & a) const;

	void	writeToStream(ostream &s) const;
	
	static String	getNameFromCharge(INT charge);
	static double	getChargeFromName(String name);


	friend ostream &	operator << (ostream &s, const Atom &);

private:
	void	assignName();
	
	
INT	nr;
double	charge;
};


inline
INT	Atom::operator == (const Atom & a) const
{	return /*a.Point::compare(b) && */ TypeClass::operator == (a);	}

inline
INT	Atom::operator != (const Atom & a) const
{	return /*a.Point::compare(b) && */ TypeClass::operator != (a);	}

#endif
