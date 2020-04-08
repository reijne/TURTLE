//***********************************************************************
//
//	Name:			Molecule.cc
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

#include "Molecule.h"

#include <iostream>
#include <stdlib.h>

#include "../../Container/GrepAwk.h"
#include "InternalCoord.h"

using namespace std;

Molecule::Molecule()
{
	atoms = NULL;
	nAtoms = 0;
}

Molecule::Molecule(istream &s)
{
	s >> nAtoms;
	atoms = new Atom * [nAtoms];
	
	for ( INT i=0 ; i<nAtoms ; i++ )
		atoms[i] = new Atom(s);
}

Molecule::Molecule(istream &s, QCProgOutputFormat::QCPackage format)
{
GrepAwk	ga(s);
	read(ga, format);
}

Molecule::Molecule(GrepAwk &ga, QCProgOutputFormat::QCPackage format)
{
	read(ga, format);
}

void	Molecule::read(GrepAwk &ga, QCProgOutputFormat::QCPackage format)
{
	switch ( format ) {
	case QCProgOutputFormat::MOLCAS:	
		ga.head();
		if ( !ga.grep("**** Cartesian Coordinates / Bohr, Angstrom ****") )
		{
			//	no atomic coordinates found
			//	return 1;
		}
		ga += 4;
		{
		INT	i = 0;
			while ( ga.getLineLength()>10 )
			{
				i++;
				ga++;
			}
			allocate(i);
			ga -= i;
		}
		for ( INT i=0 ; i<getNumberOfAtoms() ; i++ )
		{
		char	Nr[10], Label[10];
		double	x1, y1, z1, x, y, z;

			sscanf(ga.getLine().chars(), "%s %s %lf %lf %lf %lf %lf %lf", 
				Nr, Label, &x1, &y1, &z1, &x, &y, &z);

		Atom	*atom = new Atom(Label, 0, x, y, z);
			atoms[i] = atom;
			ga++;
		}
		break;

	case QCProgOutputFormat::MOLPRO:	
		ga.head();
		if ( !ga.grep("ATOMIC COORDINATES") )
		{
			//	no atomic coordinates found
			//	return 1;
		}
		ga += 4;
		{
		INT	i = 0;
			while ( ga.getLineLength()>10 )
			{
				i++;
				ga++;
			}
			allocate(i);
			ga -= i;
		}
		for ( INT i=0 ; i<getNumberOfAtoms() ; i++ )
		{
		char	Nr[10], Label[10];
		double	charge, x, y, z;

			sscanf(ga.getLine().chars(), "%s %s %lf %lf %lf %lf", 
				Nr, Label, &charge, &x, &y, &z);

		Atom	*atom = new Atom(Label, charge, 
					x*BohrToAngstrom, y*BohrToAngstrom, z*BohrToAngstrom);
			atoms[i] = atom;
			ga++;
		}
		break;

	case QCProgOutputFormat::GaussianWFN:
		ga.head();
		if ( !ga.grep("NUCLEI") )
		{
			//	no atomic coordinates found
			//	return 1;
		}
		allocate(atoi(ga.getWord(7).chars()));
		ga++;
		for ( INT i=0 ; i<getNumberOfAtoms() ; i++ )
		{
		char	Nr[10], Label[10];
		double	charge, x, y, z;
		char	dummy1[100], dummy2[100], dummy3[100], dummy4[100];

			sscanf(ga.getLine().chars(), "%s %s %s %s %lf %lf %lf %s %s %lf", 
				Label, Nr, dummy1, dummy2, &x, &y, &z, dummy3, dummy4, &charge);

		Atom	*atom = new Atom(Label, charge, 
					x*BohrToAngstrom, y*BohrToAngstrom, z*BohrToAngstrom);
//		Atom	*atom = new Atom(Label, charge, 
//					x, y, z);
			atoms[i] = atom;
			ga++;
		}
		break;

	case QCProgOutputFormat::Hondo:
		break;
	
	case QCProgOutputFormat::IntCoord:
		{
		InternalCoord	intCoord(ga);
			atoms = NULL;
			nAtoms = 0;
			*this = intCoord;
		}
		break;
	
	case QCProgOutputFormat::Gaussian:
		{
			ga.head();
			if ( !ga.grep("Input orientation:") )
				return;
				
			if ( !ga.grep("----------------------------------------------------------") )
				return;

			if ( !ga.grep("Center     Atomic              Coordinates (Angstroms)") )
				return;

			if ( !ga.grep("Number     Number             X           Y           Z") )
				return;

			ga += 2;
		Pix	line = ga.getIndex();
		INT	n = 0;
			while ( ga.getWord(1).length()<10 )
			{
				n++;
				ga++;
			}
			
			allocate(n);
			ga.setIndex(line);
			for ( INT i=0 ; i<getNumberOfAtoms() ; i++ )
			{
			char	Nr[10];
			double	x, y, z;
			INT charge;
				sscanf(ga.getLine().chars(), "%s %d %lf %lf %lf", 
					Nr, &charge, &x, &y, &z);

			String	Label ;
			Atom	*atom = new Atom(Atom::getNameFromCharge(charge)+Nr, charge, 
						x, y, z);
				atoms[i] = atom;
				ga++;
			}
		}
		break;
	
	case QCProgOutputFormat::undefinedPackage:
		break;
	}
}

Molecule::Molecule(const Molecule &m)
{
	allocate(m.getNumberOfAtoms());
	for ( INT i=0 ; i<getNumberOfAtoms() ; i++ )
	{
	Atom	*atom = new Atom(*m.atoms[i]);
		atoms[i] = atom;
	}
}


Molecule & Molecule::operator = (const Molecule &m)
{
	if ( atoms )
	{
		for ( INT i=0 ; i<nAtoms ; i++ )
			if ( atoms[i] )
				delete atoms[i];
		delete atoms;
	}
	allocate(m.getNumberOfAtoms());
	for ( INT i=0 ; i<getNumberOfAtoms() ; i++ )
	{
	Atom	*atom = new Atom(*m.atoms[i]);
		atoms[i] = atom;
	}
	return *this;
}

Molecule::Molecule(INT _nAtoms)
{
	allocate(_nAtoms);
}


void	Molecule::allocate(INT n)
{
	nAtoms = n;
	if ( nAtoms )
	{
		atoms = new Atom * [nAtoms];
		for ( INT i=0 ; i<nAtoms ; i++ )
			atoms[i] = NULL;
	}
}


Molecule::~Molecule()
{
	if ( atoms )
	{
		for ( INT i=0 ; i<nAtoms ; i++ )
			if ( atoms[i] )
				delete atoms[i];
		delete atoms;
	}
}

void	Molecule::writeToStream(ostream &s) const
{
	s << nAtoms << endl;
	for ( INT i=0 ; i<nAtoms ; i++ )
	{
		atoms[i]->writeToStream(s);
		s << endl;
	}
}

void	Molecule::getMinMaxPoint(Point &min, Point &max) const
{
	for ( INT i=0 ; i<nAtoms ; i++ )
	{
		if ( atoms[i]->getX()>max.getX() )
			max.getX() = atoms[i]->getX();
		if ( atoms[i]->getY()>max.getY() )
			max.getY() = atoms[i]->getY();
		if ( atoms[i]->getZ()>max.getZ() )
			max.getZ() = atoms[i]->getZ();
		if ( atoms[i]->getX()<min.getX() )
			min.getX() = atoms[i]->getX();
		if ( atoms[i]->getY()<min.getY() )
			min.getY() = atoms[i]->getY();
		if ( atoms[i]->getZ()<min.getZ() )
			min.getZ() = atoms[i]->getZ();
	}
}

Point	Molecule::getCenter() const
{
Point	min, max;
	getMinMaxPoint(min, max);
	return 0.5*(min+max);
}

double	Molecule::getSize(Point Center) const
{
double	max = 0;
	for ( INT i=0 ; i<nAtoms ; i++ )
	{
	double	r = (*atoms[i]-Center).getNorm2Sqr();
		if ( r>max )
			max = r;
	}
	return sqrt(max);
}


INT	operator == (const Molecule & a, const Molecule & b)
{
	if ( a.getNumberOfAtoms() != b.getNumberOfAtoms() )
		return 0;
INT	j;
	for ( INT i=0 ; i<a.getNumberOfAtoms() ; i++ )
	{	for ( j=0 ; j<b.getNumberOfAtoms() ; j++ )
		{	if ( a[i]==b[j] )
				break;
		}
		if ( j>=b.getNumberOfAtoms() )
			return 0;
	}
	return 1;
}


ostream &	operator << (ostream &s, const Molecule & a)
{
	for ( INT i=0 ; i<a.getNumberOfAtoms() ; i++ )
		s << a[i] << endl;

	return s;
}

const double Molecule::BohrToAngstrom = 0.52917641;
