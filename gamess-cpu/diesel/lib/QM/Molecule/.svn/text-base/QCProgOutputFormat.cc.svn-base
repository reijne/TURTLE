//***********************************************************************
//
//	Name:			QCProgOutputFormat.cc
//
//	Description:	
//
//	Author:			Michael Hanrath
//
//	Version:		0.0
//
//	Date:			26. Feb 1998
//
//***********************************************************************

#include "QCProgOutputFormat.h"

#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>


#include "../../Container/GrepAwk.h"

using namespace std;

QCProgOutputFormat::QCProgOutputFormat()
{
	package = undefinedPackage;
	part.part = undefinedMOLCASPart;
}


QCProgOutputFormat::QCProgOutputFormat(istream &is)
{
GrepAwk	ga(is, 1000);

	package = undefinedPackage;
	part.part = undefinedMOLCASPart;

	ga.head();
	if ( ga.grep("Gaussian, Inc.") )
	{
		if ( ga.grep("Carnegie Office Park") )
		{
			if ( ga.grep("Cite this work as") )
			{
				package = Gaussian;
				ga++;
				version = ga.getWord(2) + " " + ga.getWord(3) + " " + ga.getWord(4);
				return;
			}
		}
	}
	
	ga.head();
	if ( ga.grep("MOLCAS") )
	{
	const INT	n = 3;
	char	*partString[n] = { "R A S S C F", "HF-SCF", "SEWARD"};
	INT	format[n] = { 3, 2, 1 };
		for ( INT i=0 ; i<3 ; i++ )
		{
			ga.head();
			if ( ga.grep(partString[i]) )
			{
				package = MOLCAS;
				version = "?";
				if ( format[i]>1 )
				{
					ga.head();
					if ( ga.grep("MOLCAS version") )
						version = ga.getWord(4);
				}
				ga.head();
				if ( ga.grep("Molcas version") )
						version = ga.getWord(4);
				part.part = (MOLCASPart) format[i];
				break;
			}
		}
		return;
	}

	ga.head();
	if ( ga.grep("***  PROGRAM SYSTEM MOLPRO  ***") )
	{
		if ( ga.grep("University of Sussex") )
		{
			package = MOLPRO;
			ga.grep("Version");
			version = ga.getWord(2);
			return;
		}
	}
	
	ga.head();
	if ( ga.grep("GAUSSIAN") )
	{
		if ( ga.getWord(3)=="MOL" && ga.getWord(4)=="ORBITALS" )
		{
			if ( ga.grep("CENTRE ASSIGNMENTS") )
			{
				package = GaussianWFN;
				version = "";
				return;
			}
		}
	}


	ga.head();
	while ( !ga.illegal() )
	{
		if ( ga.getNumberOfWords()>=10 && ga.getNumberOfWords()<=11 && 
				ga.getWord(1).length()>=1 && ga.getWord(1).length()<=2 && 
				(ga.getWord(3)=="0" || ga.getWord(3)=="1") && 
				(ga.getWord(5)=="0" || ga.getWord(5)=="1") && 
				(ga.getWord(7)=="0" || ga.getWord(7)=="1") && 
				isInteger(ga.getWord(8)) &&
				isInteger(ga.getWord(9)) &&
				isInteger(ga.getWord(10)) &&
				isFloat(ga.getWord(2)) &&
				isFloat(ga.getWord(4)) &&
				isFloat(ga.getWord(6)) )
		{
			package = IntCoord;
			version = "";
			return;
		}
		ga++;
	}
}


QCProgOutputFormat::~QCProgOutputFormat()
{
}


INT	QCProgOutputFormat::isInteger(const String &s) const
{
	for ( INT i=0 ; i<(INT) s.length() ; i++ )
		if ( !isdigit(s[i]) )
			return 0;
	return 1;
}

INT	QCProgOutputFormat::isFloat(const String &s) const
{
char	*p;
	strtod(s.chars(), &p);
	return s.chars()!=p;
}


String	QCProgOutputFormat::getIdentifier() const
{
char	s[100];

	switch ( package ) {
	case undefinedPackage:
		sprintf(s, "unknown");
		break;
		
	case Hondo:
		sprintf(s, "Hondo");
		break;
		
	case MOLPRO:
		sprintf(s, "MOLPRO, Version %s", version.chars());
		break;
		
	case Gaussian:
		sprintf(s, "Gaussian, Version %s", version.chars());
		break;
		
	case GaussianWFN:
		sprintf(s, "Gaussian WFN");
		break;
		
	case MOLCAS:
	char	subformat[100];
		switch ( part.part ) {
		case undefinedMOLCASPart:
			break;
			
		case seward:
			sprintf(subformat, "SEWARD");
			break;
			
		case scf:
			sprintf(subformat, "SCF");
			break;
			
		case rasscf:
			sprintf(subformat, "RASSCF");
			break;
		}
		sprintf(s, "MOLCAS, Version %s, %s", version.chars(), subformat);
		break;
	
	case IntCoord:
		sprintf(s, "internal coordinates");
		break;
	
	case TURBOMOLE:
		sprintf(s, "TURBOMOLE, Version %s", version.chars());
		break;
	}

	return String(s);
}


QCProgOutputFormat::QCPackage	QCProgOutputFormat::getPackage() const
{	return	package;	}

QCProgOutputFormat::TPart	QCProgOutputFormat::getPart() const
{	return	part;	}


