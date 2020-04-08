//***********************************************************************
//
//	Name:			RootEnergies.cc
//
//	Description:	
//
//	Author:			Michael Hanrath
//
//	Version:		0.0
//
//	Date:			23.04.1997
//
//
//
//
//***********************************************************************

#include "RootEnergies.h"


#include <iomanip>

using namespace std;

INT	RootEnergies::storePTEnergy = 0;
INT	RootEnergies::storePTCoef = 0;


RootEnergies::RootEnergies(istream &s)
{
	e = NULL;
	c = NULL;
	s >> nRoots;
		
	if ( nRoots )
	{
		if ( storePTEnergy )
		{
			e = new EnergyType[nRoots];
			for ( INT i=0 ; i<nRoots ; i++ )
				s >> e[i];
		}
		else
		{
		EnergyType j;
			for ( INT i=0 ; i<nRoots ; i++ )
				s >> j;
		}
	}

	s >> nCSFs;
		
	if ( nCSFs )
	{
		if ( storePTCoef )
		{
			c = new PTCIType[nRoots*nCSFs];
			for ( INT i=0 ; i<nRoots*nCSFs ; i++ )
				s >> c[i];
		}
		else
		{
		EnergyType j;
			for ( INT i=0 ; i<nRoots*nCSFs ; i++ )
				s >> j;
		}
	}
}


void	RootEnergies::writeToStream(ostream &s) const
{
	if ( storePTEnergy )
	{
		s << nRoots << endl;
		for ( INT i=0 ; i<nRoots ; i++ )
			s << setw(16) << setiosflags(ios::scientific) << e[i] << " ";
		if ( nRoots )
			s << endl;
	}
	else
		s << 0 << endl;

	if ( storePTCoef )
	{
		s << nCSFs << endl;
		for ( INT i=0 ; i<nRoots*nCSFs ; i++ )
			s << setw(16) << setiosflags(ios::scientific) << c[i] << " ";
		if ( nCSFs )
			s << endl;
	}
	else
		s << 0 << endl;
}


ostream& operator<<(ostream & s, const RootEnergies &re)
{
	if ( !re.storePTEnergy && !re.storePTCoef )
		return s;

	s << "[";
	for ( INT i=0 ; i<re.nRoots ; i++ )
	{
		if ( re.storePTEnergy && re.e )
			s << setw(16) << re.e[i];
		if ( re.storePTCoef && re.c )
		{
			s << " (";
			for ( INT j=0 ; j<re.nCSFs ; j++ )
			{
				s << setw(16) << re.c[i*re.nCSFs + j];
				if ( i<re.nCSFs-1 )
					s << " ";
			}
			s << ")";
		}
		if ( i<re.nRoots-1 )
			s << ", ";
	}
	s << "]";
	return s;
}
