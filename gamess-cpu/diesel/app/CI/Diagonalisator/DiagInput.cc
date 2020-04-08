//***********************************************************************
//
//	Name:			DiagInput.cc
//
//	Description:	stores configuration input for 
//					individually selecting MR-CI
//
//	Author:			Michael Hanrath
//
//	Version:		0.0
//
//	Date:			13.09.1997
//
//
//
//
//
//***********************************************************************

#include "DiagInput.h"

#include <stdlib.h>
#include <stdio.h>
#include <sstream>

#include "../../../lib/QM/IO/Verbosity.h"
#include "../../../lib/Container/TempDir.h"

typedef double VectorType;

const VectorType DiagInput::refThreshDefault = 0.004;
const VectorType DiagInput::PTrefThreshDefault = 0.004;
const EnergyType DiagInput::ConvergenceEnergyChangeDefault = 1e-5;
const String DiagInput::ConfTreeFileNameDefault = "ConfTree.dat";
const VectorType DiagInput::ConvergenceEigenvectorChangeDefault = 1;
const INT DiagInput::MaxItersDefault = 20;
const DavidsonCI<double, double>::IterationMode  DiagInput::iterationModeDefault =
        DavidsonCI<double, double>::CI;

DiagInput::DiagInput(INT autoDetectFort31)
{
	if ( autoDetectFort31 )
	{
		MOIntegralFile = Fort31File("fort.31", Fort31RecordFormatAuto);
	Fort31File	f31(MOIntegralFile.name);
		MOIntegralFile.format = f31.format;
	}
	else
		MOIntegralFile = Fort31File("fort.31", Fort31RecordFormatUndefined);
	RefThreshold = refThreshDefault;
	PTRefThreshold = PTrefThreshDefault;
	ConvergenceEnergyChange = ConvergenceEnergyChangeDefault;
	ConvergenceEigenvectorChange = ConvergenceEigenvectorChangeDefault;
        ConfTreeFileName = ConfTreeFileNameDefault;
	rootHoming = 0;
	storePTEnergy = 0;
	storePTCoef = 0;
	MaxIters = MaxItersDefault;
	MaxStorageMem = 0;
	NumberOfRoots = 0;
	rootNumbers = NULL;
	iterationMode = iterationModeDefault;
	precision = doublePrec;
}


DiagInput::~DiagInput()
{
	if ( rootNumbers )
		free(rootNumbers);
}


DiagInput::DiagInput(const DiagInput &diag)
{
	MOIntegralFile = diag.getMOIntegralFile();
	RefThreshold = diag.getRefThreshold();
	PTRefThreshold = diag.getPTRefThreshold();
	ConvergenceEnergyChange = diag.getConvergenceEnergyChange();
	ConvergenceEigenvectorChange = diag.getConvergenceEigenvectorChange();
        ConfTreeFileName = diag.getConfTreeFileName();
	rootHoming = diag.getRootHoming();
	storePTEnergy = diag.getStorePTEnergy();
	storePTCoef = diag.getStorePTCoef();
	MaxIters = diag.getMaxIters();
	MaxStorageMem = diag.getMaxStorageMem();
	iterationMode = diag.getIterationMode();
	precision = diag.getPrecision();

	NumberOfRoots = diag.getNumberOfRoots();
	rootNumbers = (INT *) malloc(NumberOfRoots * sizeof(INT));
	for ( INT i=0 ; i<NumberOfRoots ; i++ )
		rootNumbers[i] = diag.getRootNumber(i);
}


ostream& operator<<(ostream& s, const DiagInput & diag)
{
	s << "Verbosity                    = " << verbosity << endl;
	s << "TempDir                      = " << TempDir << endl;
	s << "MOIntegralFilename           = " << diag.MOIntegralFile.name << endl;
	s << "MOIntegralFileFormat         = ";
	switch ( diag.MOIntegralFile.format )
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

	if ( diag.NumberOfRoots )
	{
		s << "Roots                        = { ";
		for ( INT i=0 ; i<diag.NumberOfRoots ; i++ )
			s << diag.getRootNumber(i) << " ";
		s << "}" << endl;
	}
	else
		s << "Roots                        = selector input" << endl;

	s << "ReferenceThreshold           = " << diag.RefThreshold << endl;
	s << "PTReferenceThreshold         = " << diag.PTRefThreshold << endl;
	s << "ConvergenceEnergyChange      = " << diag.ConvergenceEnergyChange << endl;
	s << "ConvergenceEigenvectorChange = " << diag.ConvergenceEigenvectorChange << endl;
        s << "ConfTreeFileName             = " << diag.ConfTreeFileName << endl;
	s << "MaxIters                     = " << diag.MaxIters << endl;
	s << "RootHoming                   = " << ( diag.rootHoming ? "yes" : "no") << endl;
	s << "StorePTEnergy                = " << ( diag.storePTEnergy ? "yes" : "no") << endl;
	s << "StorePTCoef                  = " << ( diag.storePTCoef ? "yes" : "no") << endl;
	s << "MaxHamiltonStorageMem        = " << (diag.MaxStorageMem >> 20) << "MB" << endl;
	s << "Precision                    = ";
	switch ( diag.precision ) {
	case DiagInput::floatPrec:
		s << "float";
		break;

	case DiagInput::doublePrec:
		s << "double";
		break;

	case DiagInput::undefinedPrec:
		s << "undefined";
		break;
	}
	s  << endl;
		
	s << "IterationMode                = ";
	switch ( diag.getIterationMode() ) {
	case DavidsonCI<double, double>::CI:
		s << "CI";
		break;
		
	case DavidsonCI<double, double>::ACPF:
		s << "ACPF";
		break;
		
	case DavidsonCI<double, double>::AQCC:
		s << "AQCC";
		break;
	}
	s  << endl;

	return s;
}


#include <FlexLexer.h>

DiagInput *LexDiag;
INT	LexDiagErrors;




std::istream& operator>>(std::istream& s, DiagInput & diag)
{
yyFlexLexer	lexer(&s);

ostream	*outStream = new std::ostringstream(std::ostringstream::in | std::ostringstream::out);


	LexDiagErrors = 0;

	diag.NumberOfRoots = 0;
	diag.rootNumbers = NULL;
	diag.RefThreshold = -1;
	diag.PTRefThreshold = -1;
	diag.ConvergenceEnergyChange = -1;
	diag.ConvergenceEigenvectorChange = -1;
        diag.ConfTreeFileName = "-1";
	diag.MaxIters = -1;
	diag.MaxStorageMem = -1;
	diag.precision = DiagInput::undefinedPrec;
	diag.iterationMode = (DavidsonCI<double, double>::IterationMode) -1;
	//diag.MOIntegralFile = Fort31File("", Fort31RecordFormatAuto);
	
	LexDiag = &diag;
	*outStream << "scanning input..." << endl << endl;
	lexer.yylex();

	if ( verbosity.isActive(Verbosity::Input) )
		outStream = &cout;
	
	if ( LexDiagErrors )
	{
		cout << endl << LexDiagErrors << " scanning errors." << endl;
		cout << "aborting." << endl;
		exit(1);
	}
	else
		*outStream << "no scanning errors." << endl;
	*outStream << endl;
	
	*outStream << "parsing input..." << endl << endl;

INT	errors = 0;

	if ( strlen(diag.MOIntegralFile.name) == 0 )
	{
		*outStream << "\"MOIntegralFilename\" is missing, defaults to \"fort.31\"." << endl;	
		diag.MOIntegralFile.name = "fort.31";
	}
		
	if ( diag.MOIntegralFile.format == Fort31RecordFormatAuto )
	{
		*outStream << "autodetecting \"MOIntegralFileFormat\"... ";	
	Fort31File	f31(diag.MOIntegralFile.name);
		diag.MOIntegralFile.format = f31.format;
		*outStream << "OK" << endl;	
		if ( f31.format==Fort31RecordFormatUndefined )
		{
			cout << "Error: unknown MOIntegralFileFormat" << endl;
			errors++;
		}
	}
	
	if ( diag.rootNumbers == NULL )
	{
		*outStream << "\"Roots\" is missing, defaults to selector input." << endl;
		diag.NumberOfRoots = 0;
		diag.rootNumbers = NULL;
	}
	
	
	if ( diag.RefThreshold == -1 )
	{
		*outStream << "\"ReferenceThreshold\" is missing, defaults to " 
			<< diag.refThreshDefault << "." << endl;
		diag.RefThreshold = diag.refThreshDefault;
	}
	
	if ( diag.PTRefThreshold == -1 )
	{
		*outStream << "\"PTReferenceThreshold\" is missing, defaults to " 
			<< diag.PTrefThreshDefault << "." << endl;
		diag.PTRefThreshold = diag.PTrefThreshDefault;
	}
	
	if ( TempDir.length()==0 )
	{
		*outStream << "\"TempDir\" is missing, defaults to \".\"." << endl;	
		TempDir = ".";
	}
	
	if ( diag.ConvergenceEnergyChange != -1 && 
			diag.ConvergenceEigenvectorChange != -1 )
	{
		cout << "Error: only one convergence criterion applicable." << endl;
		errors++;
	}
	
	if ( diag.ConvergenceEnergyChange == -1 )
	{
		*outStream << "\"ConvergenceEnergyChange\" is missing, defaults to " 
			<< diag.ConvergenceEnergyChangeDefault << "." << endl;
		diag.ConvergenceEnergyChange = diag.ConvergenceEnergyChangeDefault;
	}
	
	if ( diag.ConvergenceEigenvectorChange == -1 )
	{
		*outStream << "\"ConvergenceEigenvectorChange\" is missing, defaults to " 
			<< diag.ConvergenceEigenvectorChangeDefault << "." << endl;
		diag.ConvergenceEigenvectorChange = diag.ConvergenceEigenvectorChangeDefault;
	}

        if (diag.ConfTreeFileName == "-1" )
        {
                *outStream << "\"ConfTreeFileName\" is missing defaults to "
                        << diag.ConfTreeFileNameDefault << "." << endl;
                diag.ConfTreeFileName = diag.ConfTreeFileNameDefault;                
        }
	
	if ( diag.rootHoming == -1 )
	{
		*outStream << "\"RootHoming\" is missing, defaults to \"no\"." << endl;
		diag.rootHoming = 0;
	}

	if ( diag.storePTEnergy == -1 )
	{
		*outStream << "\"StorePTEnergy\" is missing, defaults to \"no\"." << endl;
		diag.storePTEnergy = 0;
	}

	if ( diag.storePTCoef == -1 )
	{
		*outStream << "\"StorePTCoef\" is missing, defaults to \"no\"." << endl;
		diag.storePTCoef = 0;
	}

	if ( diag.MaxIters == -1 )
	{
		*outStream << "\"MaxIters\" is missing, defaults to " 
			<< diag.MaxItersDefault << "." << endl;
		diag.MaxIters = diag.MaxItersDefault;
	}
	
	if ( diag.MaxStorageMem == -1 )
	{
		*outStream << "\"MaxStorageMem\" is missing, defaults to 0." << endl;
		diag.MaxStorageMem = 0;
	}
	
	if ( diag.iterationMode == -1 )
	{
		*outStream << "\"IterationMode\" is missing, defaults to \"CI\"." << endl;
		diag.iterationMode = diag.iterationModeDefault;
	}
	
	if ( diag.precision == DiagInput::undefinedPrec )
	{
		*outStream << "\"Precision\" is missing, defaults to \"double\"." << endl;
		diag.precision = DiagInput::doublePrec;
	}
	
	
FILE	*f;
	if (diag.MOIntegralFile.format != RIFormat)
	{
		if ( (f=fopen(diag.MOIntegralFile.name, "r"))==NULL )
		{
			cout << "MOIntegralFile.format= " << diag.MOIntegralFile.format << endl;
			cout << "Error: no file \"" << diag.MOIntegralFile.name
				<< "\"." << endl;
			errors++;
		}
		else
			fclose(f);
	}


	*outStream << endl;
	if ( errors )
	{
		cout << endl << errors << " errors.\naborting." << endl;
		exit(1);
	}
	else
	{
		*outStream << "no errors." << endl << endl;
	}
	
//	if ( !verbosity.isActive(Verbosity::Input) )
//		delete outStream;

	return s;
}

void	DiagInput::setVerbosity(const char *yytext)
{
INT	i = strlen(yytext) - 2;
	while ( i>=0 && yytext[i]!='=' )
		i--;
	i++;

        std::istringstream	s(yytext+i,std::istringstream::in | std::istringstream::out);
	verbosity = Verbosity(s);
}
