//***********************************************************************
//
//	Name:			RefSelInput.cc
//
//	Description:	stores configuration input for 
//					individually selecting MR-CI
//
//	Author:			Michael Hanrath
//
//	Version:		0.0
//
//	Date:			05.02.1997
//
//
//
//
//
//***********************************************************************

#include "RefSelInput.h"

#include <stdio.h>
#include <stdlib.h>
#include <string>

#include <sstream>

#include "../../../lib/QM/MO/MRMOs.h"
#include "../../../lib/QM/MO/MOMapping.h"

#include "../../../lib/QM/IO/Verbosity.h"

#include "../../../lib/QM/MO/Iterators/MORestriction.h"
#include "../../../lib/QM/MO/Iterators/MOEquivalence.h"
#include "../../../lib/QM/MO/MOStatistics.h"

#include "../../../lib/Container/AVLSet.h"
#include "../../../lib/QM/MO/Iterators/MOIterator.h"

using namespace std;

RefSelInput::RefSelInput()
{
	firstGuessConfs = 0;
	moEquivalence = NULL;
	mrmos = 0;
	ReferenceThreshold = 0;
	refSelMode = RefSelMode::ConfThresh;
}


RefSelInput::~RefSelInput()
{
	if ( mrmos )
		delete mrmos;
		

	if ( moEquivalence )
		delete moEquivalence;

}

RefSelInput::RefSelInput(const RefSelInput &mrconf)
{
	MOIntegralFile = mrconf.getMOIntegralFile();
	firstGuessConfs = mrconf.firstGuessConfs;
	ReferenceThreshold = mrconf.getReferenceThreshold();
	refSelMode = mrconf.getRefSelMode();

	mrmos = new MRMOs(*mrconf.getMRMOs());
	
		
	if ( mrconf.getMOEquivalence() )
		moEquivalence = new MOEquivalence(*mrconf.getMOEquivalence());
	else
		moEquivalence = NULL;
}

RefSelInput & RefSelInput::operator = (const RefSelInput &mrconf)
{
	MOIntegralFile = mrconf.getMOIntegralFile();
	firstGuessConfs = mrconf.firstGuessConfs;
	ReferenceThreshold = mrconf.getReferenceThreshold();
	refSelMode = mrconf.getRefSelMode();

	if ( mrmos )
		delete mrmos;
		

		

	mrmos = new MRMOs(*mrconf.getMRMOs());
	


	if ( moEquivalence )
		delete moEquivalence;

	if ( mrconf.getMOEquivalence() )
		moEquivalence = new MOEquivalence(*mrconf.getMOEquivalence());
	else
		moEquivalence = NULL;

	return *this;
}

#define printSet(x,s) i = x.first(); while ( i ) { s << (x)(i) << " "; (x).next(i); }

ostream& operator<<(ostream& s, const RefSelInput & mrconf)
{
	s << "Verbosity                  = " << verbosity << endl;
	s << "MOIntegralFilename         = " << mrconf.MOIntegralFile.name << endl;
	s << "MOIntegralFileFormat       = ";

	switch ( mrconf.MOIntegralFile.format )
	{
	case Fort31RecordFormatOld:
		s << "Old" << endl;
		break;
		
	case Fort31RecordFormatNew:
		s << "New" << endl;
		break;
		
	case Fort31RecordFormatTRADPT:
		s << "TRADPT" << endl;
		break;
		
	case RIFormat:
		s << "RIFormat" << endl;
		break;
		
	default:
		s << "Undefined" << endl;
		break;
	}

	if ( mrconf.moEquivalence )
		s << "MOEquivalence              = {" << mrconf.moEquivalence->getEquivalence() << "}" << endl;
	else
		s << "MOEquivalence              = none" << endl;



Pix	i = NULL;
	s << "IrReps                     = { "; printSet(mrconf.irreps,s); s << "}" << endl;


	s << "FirstGuessConfs            = " << mrconf.firstGuessConfs << endl;
	s << "ReferenceThreshold         = " << mrconf.ReferenceThreshold << endl;
	s << "RefSelMode                 = ";
	switch ( mrconf.refSelMode ) {
	case RefSelMode::ConfThresh:
		s << "ConfThresh";
		break;
	case RefSelMode::SumThresh:
		s << "SumThresh";
		break;
	}
	s << endl;

	return s;
}





#include <FlexLexer.h>

RefSelInput *LexMRConf;
INT	LexMRConfErrors;






istream& operator>>(istream& s, RefSelInput & mrconf)
{
yyFlexLexer	lexer(&s);

ostream	*outStream = new ostringstream(ostringstream::in | ostringstream::out);
	LexMRConfErrors = 0;

	mrconf.firstGuessConfs = -1;
	mrconf.ReferenceThreshold = -1;
	mrconf.refSelMode = RefSelMode::ConfThresh;
	//mrconf.MOIntegralFile = Fort31File("", Fort31RecordFormatAuto);
	
	LexMRConf = &mrconf;
	*outStream << "scanning input..." << endl << endl;
	lexer.yylex();

	if ( verbosity.isActive(Verbosity::Input) )
		outStream = &cout;

	if ( LexMRConfErrors )
	{
		cout << endl << LexMRConfErrors << " scanning errors." << endl;
		cout << "aborting." << endl;
		exit(1);
	}
	else
		*outStream << "no scanning errors." << endl;
	*outStream << endl;
	
	*outStream << "parsing input..." << endl << endl;

INT	errors = 0;



	if ( strlen(mrconf.MOIntegralFile.name) == 0 )
	{
		*outStream << "\"MOIntegralFilename\" is missing, defaults to \"fort.31\"." << endl;	
		mrconf.MOIntegralFile.name = "fort.31";
	}

	if ( !mrconf.irreps.length() )
	{
		*outStream << "\"Irreps\" is mandatory." << endl;
		errors++;
	}
	
	
	if ( mrconf.MOIntegralFile.format == Fort31RecordFormatAuto )
	{
		*outStream << "autodetecting \"MOIntegralFileFormat\"... ";	
	Fort31File	f31(mrconf.MOIntegralFile.name);
		mrconf.MOIntegralFile.format = f31.format;
		*outStream << "OK" << endl;	
		if ( f31.format==Fort31RecordFormatUndefined )
		{
			cout << "Error: unknown MOIntegralFileFormat" << endl;
			errors++;
		}
	}
	
	
	if ( mrconf.ReferenceThreshold == -1 )
	{
		*outStream << "\"ReferenceThreshold\" is missing, defaults 0.9." << endl;
		mrconf.ReferenceThreshold = 0.9;
	}
	
	
	if ( mrconf.firstGuessConfs==-1 )
	{
		*outStream << "\"FirstGuessConfs\" is missing, defaults to 1000." << endl;
		mrconf.firstGuessConfs = 1000;
	}
	
	

FILE	*f;
	if ( (f=fopen(mrconf.MOIntegralFile.name, "r"))==NULL )
	{
		cout << "Error: no file \"" << mrconf.MOIntegralFile.name
			<< "\"." << endl;
		errors++;
	}
	else
		fclose(f);



	*outStream << endl;
	if ( errors )
	{
		cout << endl << errors << " errors.\naborting." << endl;
		exit(1);
	}



	mrconf.mrmos = new MRMOs(mrconf.MOIntegralFile);

	


	*outStream << endl;









	
	*outStream << endl;
	if ( errors )
	{
		cout << endl << errors << " errors.\naborting." << endl;
		exit(1);
	}
	return s;
}




void	RefSelInput::setMOEquivalence(const char *yytext, INT auto1)
{
	if ( !yytext )
		return;
INT	i = 0;
	while ( i<(INT) strlen(yytext) && yytext[i]!='=' )
		i++;
	i++;
String	s(String(yytext+i));
	s.gsub('{',' ');
	s.gsub('}',' ');
	autoEquiv = auto1;
	if ( !auto1 )
		moEquivalence = new MOEquivalence(String(yytext+i));
}




void	RefSelInput::setVerbosity(const char *yytext)
{
INT	i = strlen(yytext) - 2;
	while ( i>=0 && yytext[i]!='=' )
		i--;
	i++;

istringstream	s(yytext+i,istringstream::in | istringstream::out);
	verbosity = Verbosity(s);
}
