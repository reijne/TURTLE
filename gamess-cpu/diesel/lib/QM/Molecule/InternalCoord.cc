//***********************************************************************
//
//	Name:			InternalCoord.cc
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

#include "InternalCoord.h"

#include <stdio.h>
#include <iostream>
#include <math.h>

#include "../../Container/GrepAwk.h"

using namespace std;


InternalCoord::InternalCoord()
{
	internal = NULL;
	n = 0;
}

InternalCoord::InternalCoord(GrepAwk &ga)
{
	read(ga);
}

InternalCoord::InternalCoord(istream &is)
{
GrepAwk	ga(is);
	read(ga);
}

void	InternalCoord::read(GrepAwk &ga)
{
	internal = NULL;
	n = 0;

	while (!ga.illegal() &&  (ga.getNumberOfWords()<10 || ga.getNumberOfWords()>11) )
		ga++;

Pix	line = ga.getIndex();
	
	n = 0;
	while ( !ga.illegal() && !(ga.getNumberOfWords()<10 || ga.getNumberOfWords()>11) )
	{
		n++;
		ga++;
	}
	ga.setIndex(line);
	
	internal = new TInternal[n];

	for ( INT i=0 ; i<n ; i++ )
	{
	char	Label[100];
		sscanf(ga.getLine().chars(),
			 "%s %lf %d %lf %d %lf %d %d %d %d",
			Label, 
			&internal[i].r, 
			&internal[i].rOpt, 
			&internal[i].phi, 
			&internal[i].phiOpt, 
			&internal[i].theta, 
			&internal[i].thetaOpt, 
			&internal[i].rel0, 
			&internal[i].rel1, 
			&internal[i].rel2
			);
		internal[i].Label = Label;
		ga++;
	}
}

InternalCoord::InternalCoord(const Molecule &molecule)
{
	n = molecule.getNumberOfAtoms();
	internal = new TInternal[n];
}

InternalCoord::operator Molecule()
{
Molecule	molecule(n);

	molecule.atoms[0] = new Atom(internal[0].Label, 
		Atom::getChargeFromName(internal[0].Label), 0, 0, 0);

	if ( n<2 )
		return molecule;
	molecule.atoms[1] = new Atom(internal[1].Label, 
		Atom::getChargeFromName(internal[1].Label), internal[1].r, 0, 0);

	if ( n<3 )
		return molecule;

double	alpha = M_PI - internal[2].phi 
		- asin(internal[2].r*sin(internal[2].phi)/internal[1].r);
	molecule.atoms[2] = new Atom(internal[2].Label, 
		Atom::getChargeFromName(internal[2].Label), 
		molecule.atoms[1]->getX()-internal[2].r*cos(alpha),
		internal[2].r*sin(alpha),
		0);

	for ( INT i=3 ; i<n ; i++ )
	{
	Point	p = getCartCoord(
			*molecule.atoms[internal[i].rel0-1],
			*molecule.atoms[internal[i].rel1-1],
			*molecule.atoms[internal[i].rel2-1],
			internal[i].r,
			internal[i].phi,
			internal[i].theta);
		molecule.atoms[i] = new Atom(internal[i].Label, 
			Atom::getChargeFromName(internal[i].Label),
			p.getX(), p.getY(), p.getZ());
	}

	return molecule;
}


InternalCoord::~InternalCoord()
{
	if ( internal )
		delete[] internal;
}

Point	InternalCoord::getCartCoord(
	Point p0, Point p1, Point p2, double r, double phi, double theta) const
{
Point	n;
	n = p1 - p0;
	n.normalize2();

double	phi1 = (180-phi)*M_PI/180.0;


Vector<double>	_center = p0 - r*cos(phi1)*n;
Point	center(_center[0], _center[1], _center[2]);
double	rk = r * sin(phi1);

Point	g;
	g = p2 - p1;

	g = g - (g*n)*n;
	g.normalize2();
	
Point	h;
	h = g % n;

	theta *= M_PI/180.0;

Point	p;
	p = center + rk*cos(theta)*g + rk*sin(theta)*h;

	return p;

}


void	InternalCoord::writeToStream(ostream &s) const
{
char	line[1000];
	for ( INT i=0 ; i<n ; i++ )
	{
		sprintf(line, "%3s %10.7f %d %14.5f %d %14.5f %d %5d %5d %5d\n",
			internal[i].Label.chars(), 
			internal[i].r, 
			internal[i].rOpt, 
			internal[i].phi, 
			internal[i].phiOpt, 
			internal[i].theta, 
			internal[i].thetaOpt, 
			internal[i].rel0, 
			internal[i].rel1, 
			internal[i].rel2
			);
		s << line;
	}
}

ostream & operator << (ostream &s, const InternalCoord &a)
{
	a.writeToStream(s);
	return s;
}
