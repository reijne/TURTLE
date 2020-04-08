//***********************************************************************
//
//	Name:			Molecule.h
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




#ifndef _MOLECULE_H
#define _MOLECULE_H

#include "../../../config.h"

#include "Atom.h"

#include "QCProgOutputFormat.h"



//FD class	ostream;
//FD class	istream;
class	GrepAwk;


class Molecule {
public:
friend	class InternalCoord;

	Molecule();
	Molecule(istream &);
	Molecule(istream &, QCProgOutputFormat::QCPackage);
	Molecule(GrepAwk &, QCProgOutputFormat::QCPackage);
	Molecule(INT _nAtoms);

	Molecule(const Molecule &);
	Molecule & operator = (const Molecule &);

	~Molecule();

	friend INT	operator == (
		const Molecule & a, const Molecule & b);
	friend INT	operator != (
		const Molecule & a, const Molecule & b)
	{	return !(a==b);	}

	INT	getNumberOfAtoms() const
	{	return nAtoms;	}
	

	Atom &	operator[](INT n) const
	{	return *atoms[n];	}	


	Point	getCenter() const;
	double	getSize(Point Center) const;
	void	getMinMaxPoint(Point &min, Point &max) const;

	friend ostream &	operator << (ostream &s, const Molecule &);
	void	writeToStream(ostream &) const;

private:
	void	read(GrepAwk &ga, QCProgOutputFormat::QCPackage format);
	void	allocate(INT n);

static const double BohrToAngstrom;

Atom	**atoms;
INT	nAtoms;
};

#endif
