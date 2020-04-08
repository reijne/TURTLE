//***********************************************************************
//
//	Name:			MRPTInput.cc
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

#include "MRPTInput.h"

#include <stdio.h>
#include <stdlib.h>
#include <string>

#include <sstream>

#include "../../../lib/QM/MO/MRMOs.h"
#include "../../../lib/QM/MO/MOMapping.h"

#include "../../../lib/QM/IO/Verbosity.h"
#include "../../../lib/Container/TempDir.h"

using namespace std;

MRPTInput::MRPTInput()
{
	nSelectionThresholds = 0;
	SelectionThreshold = NULL;
	SelectionThresholdString = NULL;
	NumberOfRoots = 0;
	rootNumbers = NULL;
	mrmos = 0;
	projectionMode = MRMPH0Matrix<double, double>::P_no0;
	inhomogenityThreshold = "";
	calcMP3 = 0;
}


MRPTInput::~MRPTInput()
{
	if ( mrmos )
		delete mrmos;
		
	if ( rootNumbers )
		free(rootNumbers);

	if ( SelectionThreshold )
		free(SelectionThreshold);

	if ( SelectionThresholdString )
	{
		for ( INT i=0 ; i<nSelectionThresholds ; i++ )
			delete SelectionThresholdString[i];
		free(SelectionThresholdString);
	}

}

MRPTInput::MRPTInput(const MRPTInput &mrconf)
{
	MOIntegralFile = mrconf.getMOIntegralFile();


	nSelectionThresholds = mrconf.getNumberOfSelectionThresholds();
	SelectionThreshold = (EnergyType *) malloc(sizeof(EnergyType) * nSelectionThresholds);
	memcpy(SelectionThreshold, mrconf.getSelectionThresholdP(),
		nSelectionThresholds*sizeof(EnergyType));

	SelectionThresholdString = (char **) malloc(sizeof(char *) * nSelectionThresholds);
	for ( INT i=0 ; i<nSelectionThresholds ; i++ )
	{
		SelectionThresholdString[i] = new char[
			strlen(mrconf.getSelectionThresholdString(i))+1];
		strcpy(SelectionThresholdString[i], mrconf.getSelectionThresholdString(i));
	}

	NumberOfRoots = mrconf.getNumberOfRoots();
	rootNumbers = (INT *) malloc(NumberOfRoots * sizeof(INT));
	for ( INT i=0 ; i<NumberOfRoots ; i++ )
		rootNumbers[i] = mrconf.getRootNumber(i);
	mrmos = new MRMOs(*mrconf.getMRMOs());
	projectionMode = mrconf.projectionMode;
	inhomogenityThreshold = mrconf.inhomogenityThreshold;
	calcMP3 = mrconf.calcMP3;
}

MRPTInput & MRPTInput::operator = (const MRPTInput &mrconf)
{
	MOIntegralFile = mrconf.getMOIntegralFile();

	if ( mrmos )
		delete mrmos;
		
	if ( rootNumbers )
		free(rootNumbers);

		
	if ( SelectionThreshold )
		free(SelectionThreshold);

	if ( SelectionThresholdString )
	{
		for ( INT i=0 ; i<nSelectionThresholds ; i++ )
			delete SelectionThresholdString[i];
		free(SelectionThresholdString);
	}



	nSelectionThresholds = mrconf.getNumberOfSelectionThresholds();
	SelectionThreshold = (EnergyType *) malloc(sizeof(EnergyType) * nSelectionThresholds);
	memcpy(SelectionThreshold, mrconf.getSelectionThresholdP(),
		nSelectionThresholds*sizeof(EnergyType));

	SelectionThresholdString = (char **) malloc(sizeof(char *) * nSelectionThresholds);
	for ( INT i=0 ; i<nSelectionThresholds ; i++ )
	{
		SelectionThresholdString[i] = new char[
			strlen(mrconf.getSelectionThresholdString(i))+1];
		strcpy(SelectionThresholdString[i], mrconf.getSelectionThresholdString(i));
	}

	NumberOfRoots = mrconf.getNumberOfRoots();
	rootNumbers = (INT *) malloc(NumberOfRoots * sizeof(INT));
	for ( INT i=0 ; i<NumberOfRoots ; i++ )
		rootNumbers[i] = mrconf.getRootNumber(i);
	mrmos = new MRMOs(*mrconf.getMRMOs());

	projectionMode = mrconf.projectionMode;
	inhomogenityThreshold = mrconf.inhomogenityThreshold;
	
	calcMP3 = mrconf.calcMP3;
	
	return *this;
}


ostream& operator<<(ostream& s, const MRPTInput & mrconf)
{
	s << "Verbosity            = " << verbosity << endl;
	s << "TempDir              = " << TempDir << endl;
	s << "MOIntegralFilename   = " << mrconf.MOIntegralFile.name << endl;
	s << "MOIntegralFileFormat = ";
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
		
	default:
		s << "Undefined" << endl;
		break;
	}


	s << "Roots                = { ";
	for ( INT i=0 ; i<mrconf.NumberOfRoots ; i++ )
		s << mrconf.getRootNumber(i) << " ";
	s << "}" << endl;


	s << "SelectionThresholds  = { ";
	for ( INT i=0 ; i<mrconf.getNumberOfSelectionThresholds() ; i++ )
		cout << mrconf.getSelectionThreshold(i) << " ";
	s << "}" << endl;

	s << "ProjectionMode       = ";
	switch ( mrconf.getProjectionMode() ) {
	case MRMPH0Matrix<double, double>::P_no0:
		s << "no0";
		break;

	case MRMPH0Matrix<double, double>::P_0Complement:
		s << "0Complement";
		break;
	}
	s << endl;

	s << "InhomogenityThreshold= " << (mrconf.getInhomogenityThreshold().length() ?
	mrconf.getInhomogenityThreshold() : String("Psi0")) << endl;

	s << "calcMP3              = " << (mrconf.getCalcMP3() ? "yes" : "no") << endl;
	return s;
}


#include <FlexLexer.h>

MRPTInput *LexMRPT;
INT	LexMRPTErrors;



static int compThresh(const void *_p1, const void *_p2)
{
const EnergyType *p1 = (const EnergyType *) _p1;
const EnergyType *p2 = (const EnergyType *) _p2;
	if ( *p1<*p2 )
		return 1;
	else
	if ( *p1>*p2 )
		return -1;
	return 0;
}


istream& operator>>(istream& s, MRPTInput & mrconf)
{
yyFlexLexer	lexer(&s);

ostream	*outStream = new ostringstream(ostringstream::in | ostringstream::out);
	LexMRPTErrors = 0;

	mrconf.NumberOfRoots = 0;
	mrconf.rootNumbers = NULL;
	mrconf.nSelectionThresholds = 0;
	mrconf.SelectionThreshold = NULL;
	mrconf.SelectionThresholdString = NULL;
	mrconf.MOIntegralFile = Fort31File("", Fort31RecordFormatAuto);
	mrconf.inhomogenityThreshold = "";
	
	LexMRPT = &mrconf;
	*outStream << "scanning input..." << endl << endl;
	lexer.yylex();

	if ( verbosity.isActive(Verbosity::Input) )
		outStream = &cout;

	if ( LexMRPTErrors )
	{
		cout << endl << LexMRPTErrors << " scanning errors." << endl;
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
	
	
	if ( TempDir.length()==0 )
	{
		*outStream << "\"TempDir\" is missing, defaults to \".\"." << endl;	
		TempDir = ".";
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
	
	if ( mrconf.rootNumbers == NULL )
	{
		*outStream << "\"Roots\" is missing, defaults to first root." << endl;
		mrconf.NumberOfRoots = 1;
		mrconf.rootNumbers = (INT *) malloc(sizeof(INT *));
		mrconf.rootNumbers[0] = 1;
	}
	
	
	if ( mrconf.nSelectionThresholds == 0 )
	{
		cout << "Error: non optional parameter \"SelectionThresholds\" missing." << endl;
		errors++;
	}

	

	if ( mrconf.inhomogenityThreshold.length() == 0 )
		*outStream << "\"InhomogenityThreshold\" is missing, defaults to Psi0." << endl;

	if ( upcase(mrconf.inhomogenityThreshold) == "PSI0" )
		mrconf.inhomogenityThreshold = "";

	

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



	qsort(mrconf.SelectionThreshold, 
		mrconf.nSelectionThresholds, sizeof(EnergyType),
		compThresh);


	mrconf.mrmos = new MRMOs(mrconf.MOIntegralFile);

	

	*outStream << endl;




	*outStream << endl;
	if ( errors )
	{
		cout << endl << errors << " errors.\naborting." << endl;
		exit(1);
	}
	else
	{
		if ( verbosity.isActive(Verbosity::MOs) )
			cout << moMapping << endl;

		mrconf.mrmos->initIntExt();
	
		if ( verbosity.isActive(Verbosity::MOs) )
		{
			cout << *mrconf.mrmos << endl;
			cout << "no errors." << endl << endl;
		}
	}
	return s;
}


void	MRPTInput::setVerbosity(const char *yytext)
{
INT	i = strlen(yytext) - 2;
	while ( i>=0 && yytext[i]!='=' )
		i--;
	i++;

istringstream	s(yytext+i,istringstream::in | istringstream::out);
	verbosity = Verbosity(s);
}
