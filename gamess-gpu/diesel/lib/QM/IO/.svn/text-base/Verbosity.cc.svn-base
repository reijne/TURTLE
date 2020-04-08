//***********************************************************************
//
//	Name:			Verbosity.cc
//
//	Description:	levels of verbosity for output
//
//	Author:			Michael Hanrath
//
//	Version:		0.0
//
//	Date:			24.02.1998
//
//
//
//
//
//***********************************************************************

#include "Verbosity.h"

#include <iostream>

#include "../../Container/String.h"

using namespace std;

Verbosity verbosity;

char	*Verbosity::SectionName[] =
	{
			"Input",
			"Integrals",
			"MOs",
			"RefGuess",
			"SGA",
			"RefMat",
			"RefMatEigenValues",
			"RefMatEigenVectors",
			"IterationBlocks",
			"WaveFunction",
			"SelectionPerRoot",
			"CacheStatistics",
			"DegenGuess",
			"DiagHist"
	};


Verbosity::Verbosity(Level level)
{
	switch ( level ) {
	case Quiet:
		flags = 0;
		break;

	case Default:
		flags = 8 | 64 | 256 | 512 | 4096;
		break;

	case Verbose:
		flags = 0xffffffff;
		break;
	}
}


Verbosity::Verbosity(INT _flags)
{
	flags = _flags;
}


Verbosity::Verbosity(istream &s)
{
	flags = 0;
String	st;
	while ( s >> st )
	{
		st = upcase(st);
		for ( INT i=0 ; i<nFlags ; i++ )
		{
		String	keyWord(Verbosity::SectionName[i]);
			keyWord = upcase(keyWord);
			if ( st==keyWord )
				setActive((Verbosity::Section) i);
		}
	}
}



ostream & operator << (ostream & s, const Verbosity & v)
{
	s << "{ ";
	for ( INT i=0 ; i<v.nFlags ; i++ )
		if ( v.isActive((Verbosity::Section) i) )
			s << Verbosity::SectionName[i] << " ";
	s << "}";
	return s;
}
