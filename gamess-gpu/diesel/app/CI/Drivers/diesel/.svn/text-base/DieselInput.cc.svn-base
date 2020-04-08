//***********************************************************************
//
//	Name:			dieselInput.cc
//
//	Description:	stores input for diesel driver 
//
//	Author:			Michael Hanrath
//
//	Version:		0.0
//
//	Date:			09.11.1998
//
//
//
//
//
//***********************************************************************

#include "DieselInput.h"

#include <stdio.h>
#include <stdlib.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <sys/wait.h>
#include <dirent.h>
#include <fcntl.h>
#include <unistd.h>

#include <iomanip>

#include <string>
#include <sstream>
#include <fstream>

#include "../../../../lib/QM/MO/MRMOs.h"
#include "../../../../lib/QM/MO/MOMapping.h"

#include "../../../../lib/QM/Configuration/ConfigurationSet.h"

#include "../../../../lib/Container/GrepAwk.h"

#include "../../../../lib/QM/IO/Verbosity.h"

#include "../../../../lib/QM/MO/Iterators/MORestriction.h"
#include "../../../../lib/QM/MO/Iterators/MOEquivalence.h"

using namespace std;

char	f31Format(Fort31RecordFormat f)
{
	switch ( f )
	{
	case Fort31RecordFormatOld:
		return 'o';
		
	case Fort31RecordFormatNew:
		return 'n';
		
	case Fort31RecordFormatTRADPT:
		return 't';

	case GAMESSC1:
		return 'g';

	case RIFormat:
		return 'r';

	default:
		cerr << "wrong f31 format" << endl;
		exit(1);
	}
}


DieselInput::DieselInput()
{
	MOIntegralFile = Fort31File("fort.31", Fort31RecordFormatUndefined);
	mrmos = 0;
	moRestriction = NULL;
	moEquivalence = NULL;
	MaxHamiltonStorageMem = 0;

	NumberOfElectrons = 0;
	ExcitationLevel = 2;
	selectInternal = 0;
	autoRef = 0;
	autoEquiv = 0;
	estimationMode = EnergyMap::EpsteinNesbet;
	rootHoming = 0;
	usePreviousSelection = 0;
	homogeneousSelection = 0;
	
	RefThreshold = 0.004;
	PTRefThreshold = 0.004;
	NatOrbRefThreshold = 0.004;
	ConvergenceEnergyChange = 1e-5;
	ConvergenceEigenvectorChange = 1;
	useNaturalOrbitals = 0;
	averagedNaturalOrbitals = 0;

//	projectionMode = ;

	firstGuessConfs = 1000;

	activeSpaceExcitationLevel = 2;
	maxRefOpenShells = 4;
	activeReferenceThreshold = 0;
	refSelMode = RefSelMode::ConfThresh;

	nProc = 0;
	MaxIters = 20;
	maxRefGenIters = 6;
}


DieselInput::~DieselInput()
{
	if ( mrmos )
		delete mrmos;

	if ( moRestriction )
		delete moRestriction;

	if ( moEquivalence )
		delete moEquivalence;
}

DieselInput::DieselInput(const DieselInput &mrconf)
{
	MOIntegralFile = mrconf.getMOIntegralFile();

	roots = mrconf.roots;
	preScanRoots = mrconf.preScanRoots;
	refSelMode = mrconf.refSelMode;

	mrmos = new MRMOs(*mrconf.getMRMOs());
	if ( mrconf.getMORestriction() )
		moRestriction = new MORestriction(*mrconf.getMORestriction());
	else
		moRestriction = NULL;

	if ( mrconf.getMOEquivalence() )
		moEquivalence = new MOEquivalence(*mrconf.getMOEquivalence());
	else
		moEquivalence = NULL;
}

DieselInput & DieselInput::operator = (const DieselInput &mrconf)
{
	MOIntegralFile = mrconf.getMOIntegralFile();
	refSelMode = mrconf.refSelMode;
	roots = mrconf.roots;
	preScanRoots = mrconf.preScanRoots;

	if ( mrmos )
		delete mrmos;
		
	mrmos = new MRMOs(*mrconf.getMRMOs());

	if ( moRestriction )
		delete moRestriction;

	if ( moEquivalence )
		delete moEquivalence;

	if ( mrconf.getMORestriction() )
		moRestriction = new MORestriction(*mrconf.getMORestriction());
	else
		moRestriction = NULL;

	if ( mrconf.getMOEquivalence() )
		moEquivalence = new MOEquivalence(*mrconf.getMOEquivalence());
	else
		moEquivalence = NULL;
	
	return *this;
}

#define printSet(x,s) i = x.first(); while ( i ) { s << (x)(i) << " "; (x).next(i); }
#define printVector(x,s) for ( vector<INT>::const_iterator _k=(x).begin() ; _k!=(x).end() ; ++_k ) s << *_k << " ";

ostream& operator<<(ostream& s, const DieselInput & diesel)
{
	s << "Verbosity                        = " << verbosity << endl;
	s << "TempDir                          = " << diesel.TempDir << endl;
	s << "MOLCASRootDir                    = " << diesel.MOLCASRootDir << endl;
	s << "MOIntegralFilename               = " << diesel.MOIntegralFile.name << endl;
	s << "MOIntegralFileFormat             = ";
	switch ( diesel.MOIntegralFile.format )
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
		
	case GAMESSC1:
		s << "GAMESSC1" << endl;
		break;
		
	case RIFormat:
		s << "RIFormat" << endl;
		break;
		
	default:
		s << "Undefined" << endl;
		break;
	}


	if ( diesel.moRestriction )
		s << "MORestrictions                   = {" << diesel.moRestriction->getRestriction() << "}" << endl;
	else
		s << "MORestrictions                   = none" << endl;

	if ( diesel.moEquivalence )
		s << "MOEquivalence                    = {" << diesel.moEquivalence->getEquivalence() << "}" << endl;
	else
		s << "MOEquivalence                    = none" << endl;


	s << "NumberOfElectrons                = " << diesel.NumberOfElectrons << endl;
	s << "ExcitationLevel                  = " << diesel.ExcitationLevel << endl;

	s << "selectInternal                   = " << ( diesel.selectInternal ? "yes" : "no") << endl;

	s << "EstimationMode                   = ";
	switch ( diesel.estimationMode )
	{
	case EnergyMap::EpsteinNesbet:
		s << "EpsteinNesbet" << endl;
		break;
		
	case EnergyMap::Wenzel:
		s << "Wenzel" << endl;
		break;
		
	default:
		s << "Undefined" << endl;
		break;
	}

	{
	INT jj = 0;
		for ( vector<vector<INT> >::const_iterator j=diesel.roots.begin() ;
			j!=diesel.roots.end() ; ++j )
		{
			if ( j->size() )
			{
				s << "Roots" << jj << "                          = { "; 
				printVector(*j, s);
				s << "}" << endl;
			}
			jj++;
		}
	}

Pix	i = NULL;
	s << "preScanRoots                     = " << diesel.preScanRoots << endl;
	s << "SelectionThresholds              = { "; printSet(diesel.thresholds,s); s << "}" << endl;
	s << "HomogeneousSelection             = " << ( diesel.homogeneousSelection ? "yes" : "no") << endl;
	s << "IrReps                           = { "; printSet(diesel.irreps,s); s << "}" << endl;
	s << "Multiplicities                   = { "; printSet(diesel.multiplicities,s); s << "}" << endl;
	s << "selectNthExcitation              = { "; printSet(diesel.nthExcitation,s); s << "}" << endl;
	s << "RefConfs                   = " // << diesel.Refs.n
		<< endl;
	s << "{" << endl;
	((DieselInput *) &diesel)->Refs.setNumbering(1);
	s << diesel.Refs;
	((DieselInput *) &diesel)->Refs.setNumbering(0);
	s << "}" << endl;
	s << "AnnihilatorSpace                 = { "; printSet(diesel.annihilatorSpace,s); s << "}" << endl;
	s << "CreatorSpace                     = { "; printSet(diesel.creatorSpace,s); s << "}" << endl;
	s << "activeSpaceExcitationLevel       = " << diesel.activeSpaceExcitationLevel << endl;
	s << "maxRefOpenShells                 = " << diesel.maxRefOpenShells << endl;
	s << "activeReferenceThreshold         = " << diesel.activeReferenceThreshold << endl;
	s << "fullMRCIExtrapolation            = { "; printSet(diesel.fullMRCIExtrapolation,s); s << "}" << endl;
	s << "propertyTreshholds               = { "; printSet(diesel.propertyThresholds,s); s << "}" << endl;
	s << "orbitalFile                      = " << diesel.orbitalFile << endl;
	s << "store                            = { "; printSet(diesel.store,s); s << "}" << endl;
	
	s << "ReferenceThreshold               = " << diesel.RefThreshold << endl;
	s << "RefSelMode                       = ";
	switch ( diesel.refSelMode ) {
	case RefSelMode::ConfThresh:
		s << "ConfThresh";
		break;
	case RefSelMode::SumThresh:
		s << "SumThresh";
		break;
	}
	s << endl;
	s << "PTReferenceThreshold             = " << diesel.PTRefThreshold << endl;
	s << "useNaturalOrbitals               = " << ( diesel.useNaturalOrbitals ? "yes" : "no") << endl;
	s << "NatOrbReferenceThreshold         = " << diesel.NatOrbRefThreshold << endl;
	s << "averagedNaturalOrbitals          = " << ( diesel.averagedNaturalOrbitals ? "yes" : "no") << endl;
	s << "nProcessors                      = " << diesel.nProc << endl;
	s << "ConvergenceEnergyChange          = " << diesel.ConvergenceEnergyChange << endl;
	s << "ConvergenceEigenvectorChange     = " << diesel.ConvergenceEigenvectorChange << endl;
	s << "MaxDavidsonIters                 = " << diesel.MaxIters << endl;
	s << "RootHoming                       = " << ( diesel.rootHoming ? "yes" : "no") << endl;
	s << "usePreviousSelection             = " << ( diesel.usePreviousSelection ? "yes" : "no") << endl;
	s << "maxRefGenIters                   = " << diesel.maxRefGenIters << endl;
	s << "MaxHamiltonStorageMem            = " << (diesel.MaxHamiltonStorageMem >> 20) << "MB" << endl;

	s << "MRPTInhomogenityThreshold        = " << diesel.MRPTInhomogenityThreshold << endl;
	s << "MRPTSelectionThresholds          = { "; printSet(diesel.MRPTthresholds,s); s << "}" << endl;

	s << endl;

	return s;
}


#include <FlexLexer.h>

DieselInput *LexDiesel;
INT	LexDieselErrors;



istream& operator>>(istream& s, DieselInput & diesel)
{
yyFlexLexer	lexer(&s);

ostream	*outStream = new ostringstream(ostringstream::in | ostringstream::out);
	LexDieselErrors = 0;

	diesel.MOIntegralFile = Fort31File("", Fort31RecordFormatAuto);
	diesel.MaxHamiltonStorageMem = 0;
	diesel.PTRefThreshold = -1;
	diesel.NatOrbRefThreshold = -1;
	diesel.firstGuessConfs = -1;
	diesel.rootHoming = -1;
	diesel.usePreviousSelection = -1;
	diesel.homogeneousSelection = -1;
	diesel.maxRefOpenShells = -1;
	diesel.activeReferenceThreshold = -1;
	diesel.refSelMode = RefSelMode::ConfThresh;
	diesel.roots.resize(8);
	diesel.preScanRoots = -1;

	LexDiesel = &diesel;
	*outStream << "scanning input..." << endl << endl;
	lexer.yylex();

	if ( verbosity.isActive(Verbosity::Input) )
		outStream = &cout;

	if ( LexDieselErrors )
	{
		cerr << endl << LexDieselErrors << " scanning errors." << endl;
		cerr << "aborting." << endl;
		exit(1);
	}
	else
		*outStream << "no scanning errors." << endl;
	*outStream << endl;
	
	*outStream << "parsing input..." << endl << endl;

INT	errors = 0;

	if ( diesel.MOLCASRootDir.length() == 0 )
	{
		cerr << "Error: non optional parameter \"MOLCASRootDir\" missing." << endl;
		errors++;
	}
	else
		diesel.MOLCASRootDir += "/";


	if ( diesel.TempDir.length() == 0 )
	{
		*outStream << "\"TempDir\" is missing, defaults to \".\"." << endl;	
		diesel.TempDir = ".";
	}
	{
	char	num[100];
		sprintf(num, "%d", getpid());
		diesel.TempDir += String("/") + String(num);
	}
	
	
	if ( strlen(diesel.MOIntegralFile.name) == 0 )
	{
		*outStream << "\"MOIntegralFilename\" is missing, defaults to \"fort.31\"." << endl;	
		diesel.MOIntegralFile.name = diesel.MOLCASRootDir+"fort.31";
	}
	
	
	if ( diesel.MOIntegralFile.format == Fort31RecordFormatAuto )
	{
		*outStream << "autodetecting \"MOIntegralFileFormat\"... ";	
	Fort31File	f31(diesel.MOIntegralFile.name);
		diesel.MOIntegralFile.format = f31.format;
		*outStream << "OK" << endl;	
		if ( f31.format==Fort31RecordFormatUndefined )
		{
			cerr << "Error: unknown MOIntegralFileFormat" << endl;
			errors++;
		}
	}
	
	if ( !diesel.rootsDefault.size() )
	{
		*outStream << "\"RootsDefault\" is missing, defaults to first root." << endl;
		diesel.rootsDefault.push_back(1);
	}


	
	if ( diesel.preScanRoots==-1 )
	{
		*outStream << "\"preScanRoots\" is missing, defaults to 10." << endl;
		diesel.preScanRoots = 10;
	}
	
	
/*	if ( !diesel.thresholds.length() )
	{
		cerr << "Error: non optional parameter \"SelectionThresholds\" missing." << endl;
		errors++;
	}
*/

	if ( diesel.maxRefOpenShells==-1 )
	{
		*outStream << "\"MaxRefOpenShells\" is missing, defaults to 4" << endl;
		diesel.maxRefOpenShells = 4;
	}

	if ( diesel.activeReferenceThreshold==-1 )
	{
		*outStream << "\"ActiveReferenceThreshold\" is missing, defaults to 0" << endl;
		diesel.activeReferenceThreshold = 0;
	}

	
	if ( diesel.rootHoming==-1 )
	{
		*outStream << "\"RootHoming\" is missing, defaults to \"no\"." << endl;
		diesel.rootHoming = 0;
	}

	if ( diesel.usePreviousSelection==-1 )
	{
		*outStream << "\"usePreviousSelection\" is missing, defaults to \"no\"." << endl;
		diesel.usePreviousSelection = 0;
	}

	if ( diesel.homogeneousSelection==-1 )
	{
		*outStream << "\"homogeneousSelection\" is missing, defaults to \"no\"." << endl;
		diesel.homogeneousSelection = 0;
	}

	if ( diesel.firstGuessConfs==-1 )
	{
		*outStream << "\"FirstGuessConfs\" is missing, defaults to 1000." << endl;
		diesel.firstGuessConfs = 1000;
	}

	if ( diesel.PTRefThreshold == -1 )
	{
		*outStream << "\"PTRefThreshold\" is missing, defaults to \"RefThreshold\"." << endl;
		diesel.PTRefThreshold = diesel.RefThreshold;
	}

	if ( diesel.NatOrbRefThreshold == -1 )
	{
		*outStream << "\"NatOrbRefThreshold\" is missing, defaults to \"RefThreshold\"." << endl;
		diesel.NatOrbRefThreshold = diesel.RefThreshold;
	}


	if ( diesel.useNaturalOrbitals && diesel.propertyThresholds.length()>0 )
	{
		cerr << "Error: property calculation does work with natural orbitals." << endl;
		errors++;
	}
	
	

FILE	*f;
	if ( (f=fopen(diesel.MOIntegralFile.name, "r"))==NULL )
	{
		cerr << "Error: no file \"" << diesel.MOIntegralFile.name
			<< "\"." << endl;
		errors++;
	}
	else
		fclose(f);



	*outStream << endl;
	if ( errors )
	{
		cerr << endl << errors << " errors.\naborting." << endl;
		exit(1);
	}




	diesel.mrmos = new MRMOs(diesel.MOIntegralFile);
//	cout << diesel.MOIntegralFile.name << " " << diesel.MOIntegralFile.format << endl;
//	cout << *diesel.mrmos << endl;

	*outStream << endl;

	{
	Pix	iIrrep = diesel.irreps.first();
		while ( iIrrep )
		{
			if ( !diesel.roots[diesel.irreps(iIrrep)].size() )
				diesel.roots[diesel.irreps(iIrrep)] = diesel.rootsDefault;
			diesel.irreps.next(iIrrep);
		}
	}

	*outStream << endl;
	if ( errors )
	{
		cerr << endl << errors << " errors.\naborting." << endl;
		exit(1);
	}
	else
	{
		if ( verbosity.isActive(Verbosity::MOs) )
			cout << moMapping << endl;

		diesel.mrmos->initIntExt();
	
		if ( verbosity.isActive(Verbosity::MOs) )
		{
			cout << *diesel.mrmos << endl;
			cout << "no errors." << endl << endl;
		}
	}
	return s;
}


void	DieselInput::setVerbosity(const char *yytext)
{
INT	i = strlen(yytext) - 2;
	while ( i>=0 && yytext[i]!='=' )
		i--;
	i++;

istringstream	s(yytext+i,istringstream::in | istringstream::out);
	verbosity = Verbosity(s);
}



void	DieselInput::writeRefGuessInput(
	String filename,
	vector<INT>	roots,
	INT multiplicity) const
{
ofstream	s(filename);

	s << "NumberOfElectrons                = " << NumberOfElectrons << endl;
	s << "Roots                            = { "; printVector(roots,s); s << "}" << endl;
	s << "Multiplicity                     = " << multiplicity << endl;
//	s << "FirstGuessConfs                  = " << firstGuessConfs << endl;
//	s << "ReferenceThreshold               = " << RefThreshold << endl;
	s << endl;
}

void	DieselInput::writeRefSelInput(
	String filename, SLList<INT> irreps) const
{
ofstream	s(filename);


	s << "MOIntegralFileFormat             = ";
	switch ( MOIntegralFile.format )
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
		
	case GAMESSC1:
		s << "GAMESSC1" << endl;
		break;
		
	case RIFormat:
		s << "RIFormat" << endl;
		break;
		
	default:
		s << "Undefined" << endl;
		break;
	}


/*	{
	INT jj = 0;
		for ( vector<vector<INT> >::const_iterator j=roots.begin() ;
			j!=roots.end() ; ++j )
		{
			s << "Roots" << jj++ << "                          = { "; 
			printVector(*j, s);
			s << "}" << endl;
		}
	}
*/

//	s << "FirstGuessConfs                  = " << firstGuessConfs << endl;
	s << "ReferenceThreshold               = " << RefThreshold << endl;
	s << "RefSelMode                       = ";
	switch ( refSelMode ) {
	case RefSelMode::ConfThresh:
		s << "ConfThresh";
		break;
	case RefSelMode::SumThresh:
		s << "SumThresh";
		break;
	}
	s << endl;

	if ( autoEquiv )
		s << "MOEquivalence                    = auto" << endl;
	else
	{
		if ( moEquivalence )
			s << "MOEquivalence                    = {" << moEquivalence->getEquivalence() << "}" << endl;
		else
			s << "MOEquivalence                    = none" << endl;
	}
Pix	i = NULL;
	s << "IrReps                           = { "; printSet(irreps,s); s << "}" << endl;
	s << endl;
}


void	DieselInput::writeSelectorInput(
	String filename,
	vector<INT>	roots,
	SLList<INT>	nthExcitation,
	INT multiplicity,
	INT irrep,
	const SLList<String> &thresholds,
	const ConfigurationSet &refConfs,
	const ConfigurationSet &PTRefConfs,
	bool supressActive) const
{
ofstream	s(filename);
Verbosity	v(verbosity);
	v.setActive(Verbosity::RefMatEigenVectors);
	s << "Verbosity                        = " << v << endl;
	s << "MOIntegralFileFormat             = ";
	switch ( MOIntegralFile.format )
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
		
	case GAMESSC1:
		s << "GAMESSC1" << endl;
		break;
		
	case RIFormat:
		s << "RIFormat" << endl;
		break;
		
	default:
		s << "Undefined" << endl;
		break;
	}


	if ( moRestriction )
		s << "MORestrictions                   = {" << moRestriction->getRestriction() << "}" << endl;
	else
		s << "MORestrictions                   = none" << endl;


/*	if ( autoEquiv )
		s << "MOEquivalence                    = auto" << endl;
	else
	{
		if ( moEquivalence )
			s << "MOEquivalence                    = {" << moEquivalence->getEquivalence() << "}" << endl;
		else
			s << "MOEquivalence                    = none" << endl;
	}
*/

Pix	i = NULL;
	s << "NumberOfElectrons                = " << NumberOfElectrons << endl;
	if ( activeReferenceThreshold!=0 && !supressActive )
		s << "ExcitationLevel                  = 0" << endl;
	else
	{
		s << "ExcitationLevel                  = " << ExcitationLevel << endl;
		s << "selectNthExcitation              = { "; printSet(nthExcitation,s); s << "}" << endl;
	}

	s << "selectInternal                   = " << ( selectInternal ? "yes" : "no") << endl;

	s << "EstimationMode                   = ";
	switch ( estimationMode )
	{
	case EnergyMap::EpsteinNesbet:
		s << "EpsteinNesbet" << endl;
		break;
		
	case EnergyMap::Wenzel:
		s << "Wenzel" << endl;
		break;
		
	default:
		s << "Undefined" << endl;
		break;
	}

	s << "Roots                            = { "; printVector(roots,s); s << "}" << endl;
	s << "SelectionThresholds              = { "; printSet(thresholds,s); s << "}" << endl;
	s << "IrRep                            = " << irrep << endl;
	s << "Multiplicity                     = " << multiplicity << endl;

	s << "RefConfs                         = ";
	if ( refConfs.length() )
	{
		((ConfigurationSet &) refConfs).setNumbering(1);
		s << "{ " << endl << refConfs << "}" << endl;
		((ConfigurationSet &) refConfs).setNumbering(0);
	}
	else
		s << "auto" << endl;

	if ( PTRefConfs!=refConfs )
	{
		s << "PTRefConfs                         = ";
		((ConfigurationSet &) PTRefConfs).setNumbering(1);
		s << "{ " << endl << PTRefConfs << "}" << endl;
		((ConfigurationSet &) PTRefConfs).setNumbering(0);
	}

	if ( !supressActive )
	{
		s << "AnnihilatorSpace                 = { "; printSet(annihilatorSpace,s); s << "}" << endl;
		s << "CreatorSpace                     = { "; printSet(creatorSpace,s); s << "}" << endl;
		s << "activeSpaceExcitationLevel       = " << activeSpaceExcitationLevel << endl;
	}
	s << "maxRefOpenShells                 = " << maxRefOpenShells << endl;
//	s << "activeReferenceThreshold         = " << activeReferenceThreshold << endl;
	s << endl;
}





void	DieselInput::writeMRPTInput(
	String filename,
	vector<INT>	roots,
	const SLList<String> &thresholds,
	INT mp3) const
{
ofstream	s(filename);

//	s << "Verbosity                        = " << verbosity << endl;
	s << "TempDir                          = " << TempDir << endl;
	s << "MOIntegralFileFormat             = ";
	switch ( MOIntegralFile.format )
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
		
	case GAMESSC1:
		s << "GAMESSC1" << endl;
		break;
		
	case RIFormat:
		s << "RIFormat" << endl;
		break;
		
	default:
		s << "Undefined" << endl;
		break;
	}



Pix	i = NULL;
	s << "Roots                            = { "; printVector(roots,s); s << "}" << endl;
	s << "InhomogenityThreshold            = " << MRPTInhomogenityThreshold << endl;
	s << "SelectionThresholds              = { "; printSet(MRPTthresholds,s); s << "}" << endl;
	s << "calcMP3                          = " << (mp3 ? "yes" : "no") << endl;

	s << endl;
}


void	DieselInput::writeDiagonalisatorInput(
	String filename,
	VectorType	RefThreshold,
	VectorType	PTRefThreshold,
	vector<INT>	roots
	) const
{
ofstream	s(filename);

//	s << "Verbosity                        = " << verbosity << endl;
	s << "TempDir                          = " << TempDir << endl;
	s << "MOIntegralFileFormat             = ";
	switch ( MOIntegralFile.format )
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
		
	case GAMESSC1:
		s << "GAMESSC1" << endl;
		break;
		
	case RIFormat:
		s << "RIFormat" << endl;
		break;
		
	default:
		s << "Undefined" << endl;
		break;
	}


	s << "ReferenceThreshold               = " << RefThreshold << endl;
	s << "PTReferenceThreshold             = " << PTRefThreshold << endl;
	if ( ConvergenceEigenvectorChange==1 )
		s << "ConvergenceEnergyChange          = " << ConvergenceEnergyChange << endl;
	else
		s << "ConvergenceEigenvectorChange     = " << ConvergenceEigenvectorChange << endl;
	s << "MaxIters                         = " << MaxIters << endl;
	s << "Roots                            = { "; printVector(roots,s); s << "}" << endl;
	s << "RootHoming                       = " << ( rootHoming ? "yes" : "no") << endl;
	s << "MaxHamiltonStorageMem            = " << (MaxHamiltonStorageMem >> 20) << "MB" << endl;

	s << endl;

}

void	DieselInput::createReferenceSpace(INT multiplicity)
{
/*
	if ( autoRef )
	{
		cerr << "creating reference space:" << endl;

		writeRefGuessInput("refguess.in", roots, multiplicity);
		
		String	command(getenv("DIESEL_EXE_DIR"));
			command += String("/refguess <refguess.in >refguess.out");
		system(command.chars());


		cerr << "reference space generation completed" << endl;



	Pix	iIrrep = irreps.first();
		while ( iIrrep )
		{
		char	buf[100];
			cerr << "\tirrep=" << irreps(iIrrep) << endl;
			sprintf(buf, "%d", irreps(iIrrep));
			mkdir(buf, S_IRUSR | S_IRGRP | S_IROTH | S_IXUSR | S_IXGRP | S_IXOTH | S_IWUSR);
			chdir(buf);

		ConfigurationSet	refConfs, PTRefConfs;

			sprintf(buf, "../refs.%d", irreps(iIrrep));
		ifstream	selInAll(buf);
		GrepAwk	ga(selInAll);

			do {
			const INT max = 10000;
			char	buf[max];
			strstream	sbuf(buf, max);
			string	s = ga.getLine().chars();
				for ( INT i=1 ; i<=ga.getNumberOfWords() ; i++ )
				{
					if ( ga.getWord(i)=="#" )
						break;
					sbuf << ga.getWord(i) << " ";
				}
				sbuf << endl;
			Configuration<MOType>	conf(sbuf);
				refConfs.add(conf);

				ga++;
			} while ( ga.getNumberOfWords()>2 );




		FILE	*fd;
			if ( !(fd=fopen("sel.in.all", "r")) )
				writeSelectorInput("sel.in.all", 
					roots, nthExcitation, multiplicity, irreps(iIrrep), thresholds, refConfs, PTRefConfs);
			else
			{
				cerr << "\t\tusing previously generated reference space" << endl;
				fclose(fd);
			}

			irreps.next(iIrrep);
			chdir("..");
		}
	}
	else
	{
	Pix	iIrrep = irreps.first();
		while ( iIrrep )
		{
		char	buf[100];
			cerr << "\tirrep=" << irreps(iIrrep) << endl;
			sprintf(buf, "%d", irreps(iIrrep));
			mkdir(buf, S_IRUSR | S_IRGRP | S_IROTH | S_IXUSR | S_IXGRP | S_IXOTH | S_IWUSR);
			chdir(buf);

		ConfigurationSet	refConfs(Refs), PTRefConfs;
			if ( creatorSpace.length() )
				cerr << "\t\tusing active space references" << endl;
			else
			{
				cerr << "\t\tusing given references" << endl;
				refConfs.projectOnIrrep(irreps(iIrrep), *mrmos);
			}
			writeSelectorInput("sel.in.all", 
				roots, nthExcitation, multiplicity, irreps(iIrrep), thresholds, refConfs, PTRefConfs);

			irreps.next(iIrrep);
			chdir("..");
		}
	}

*/









	cerr << "\tcreating reference spaces:" << endl;


Pix	iIrrep = irreps.first();
	while ( iIrrep )
	{
	char	buf[100];
	INT	irrep = irreps(iIrrep);
		cerr << "\tirrep=" << irrep << endl;
		sprintf(buf, "%d", irreps(iIrrep));
		mkdir(buf, S_IRUSR | S_IRGRP | S_IROTH | S_IXUSR | S_IXGRP | S_IXOTH | S_IWUSR);
		chdir(buf);
		symlink(String(MOIntegralFile.name), "fort.31");

	FILE	*fd;


		if ( autoRef || activeReferenceThreshold!=0 )
		{
			if ( !(fd=fopen("sel.in.all", "r"))  )
			{
			ConfigurationSet	oldRefConfs, oldPTRefConfs, refConfs, PTRefConfs;



				for ( INT iter=0 ; iter<(activeReferenceThreshold!=0 ? 1 : maxRefGenIters) ; iter++ )
				{
				SLList<String>	genThresholds;
				String	genThreshold;
				vector<INT>	genRoots;
				SLList<INT>	genNthExcitation;
					if ( iter==0 )
					{
						cerr << "\t\tmaking first guess" << endl;
						if ( activeReferenceThreshold!=0 )
						{
							genThreshold = "1";
							refConfs = Refs;
							genRoots = roots[irrep];
						}
						else
						{
							genThreshold = "0.0001";
						INT maxRoot = 0;
							for ( vector<INT>::const_iterator i=roots[irrep].begin() ; i!=roots[irrep].end() ; ++i )
								if ( maxRoot<*i )
									maxRoot = *i;
							for ( INT j=1 ; j<=maxRoot+preScanRoots ; j++ )
								genRoots.push_back(j);
						}


						writeDiagonalisatorInput("diag.in.RefGen", 
							RefThreshold, PTRefThreshold, genRoots);

					INT r = 1;
						genNthExcitation.append(r);

					}
					else
					{
						genThreshold = "0.0001";
						cerr << "\t\titeration #" << iter << endl;
					INT doNatOrb = 0;
						genRoots = roots[irrep];
//						genRoots = (!doNatOrb && useNaturalOrbitals) ? roots : preScanRoots;
						if ( doNatOrb )
							writeDiagonalisatorInput("diag.in.RefGen", 
								NatOrbRefThreshold, NatOrbRefThreshold, genRoots);
						else
							writeDiagonalisatorInput("diag.in.RefGen", 
								RefThreshold, PTRefThreshold, genRoots);
						genNthExcitation = nthExcitation;
					}

					genThresholds.append(genThreshold);

				char	selin[100];
					sprintf(selin, "sel.in.%d", iter);
				char	gen[100];
					sprintf(gen, "genspace.%d", iter);
					writeSelectorInput(selin, 
						genRoots,
						genNthExcitation,
						multiplicity, irrep, genThresholds, refConfs, PTRefConfs);

					if ( !fork() ) //	child
					{
					String	command(getenv("DIESEL_EXE_DIR"));
						command += String("/sel <") + String(selin) + String(" >") + String(gen);
						execl("/bin/SH", "SH", "-c", command.chars(), NULL);
						cerr << "error executing " << command << endl;
						exit(1);
					}
					wait(0);

					unlink(String("ConfTree.dat"));
					symlink(String("ConfTree.dat.")+genThreshold, "ConfTree.dat");
					if ( !fork() ) //	child
					{
					String	command(getenv("DIESEL_EXE_DIR"));
						command += String("/diag -i <diag.in.RefGen >>") + String(gen);
						if ( nProc>1 )
						{
						char	b[1000];
							sprintf(b, " -p %d", nProc);
							command += b;
						}
						execl("/bin/SH", "SH", "-c", command.chars(), NULL);
						cerr << "error executing " << command << endl;
						exit(1);
					}
					wait(0);

					rename(String("ConfTree.dat.")+genThreshold, "ConfTree.dat.ref");
					rename("Eigenvectors.dat", "Eigenvectors.dat.ref");

					oldRefConfs = refConfs;
					refConfs.clear();

					oldPTRefConfs = PTRefConfs;
					PTRefConfs.clear();

				SLList<INT>	ir;
					ir.append(irreps(iIrrep));
					writeRefSelInput("refsel.in", ir);

					if ( !fork() ) //	child
					{
					String	command(getenv("DIESEL_EXE_DIR"));
					char	buf[10];
						sprintf(buf, "%d", roots[irrep].size());
						command += String("/refsel ") + buf + String(" <refsel.in >>")+gen;
						execl("/bin/SH", "SH", "-c", command.chars(), NULL);
						cerr << "error executing " << command << endl;
						exit(1);
					}
					wait(0);


				ifstream	is("refs.out");
				GrepAwk	ga(is);
					while ( !ga.illegal() )
					{
					String	buf(ga.getLine());
					istringstream	s(buf.chars(),istringstream::in | istringstream::out);
					Configuration<MOType>	conf(s);
						refConfs.add(conf);
						ga++;
					}


					if ( refConfs==oldRefConfs )
						break;
				}


//===========================================================================


			}
			else
			{
				cerr << "\t\tusing previously generated reference space" << endl;
				fclose(fd);
			}
		}
		else
		{ 
			if ( activeReferenceThreshold==0 )
			{
				if ( creatorSpace.length() )
				{
					cerr << "\t\tusing active space references" << endl;
				ConfigurationSet	PTRefConfs;
					writeSelectorInput("sel.in.all", 
						roots[irrep], nthExcitation, multiplicity, irrep, thresholds, Refs, PTRefConfs);
				}
				else
				{
					cerr << "\t\tusing given references" << endl;
				ConfigurationSet	PTRefConfs;
					writeSelectorInput("sel.in.all", 
						roots[irrep], nthExcitation, multiplicity, irrep, thresholds, Refs, PTRefConfs);
				}
			}
		}
		
		
		chdir("..");
		irreps.next(iIrrep);
	}
		

	if ( autoRef || activeReferenceThreshold!=0 )
	{
		symlink(String(MOIntegralFile.name), "fort.31");
		writeRefSelInput("refsel.in", irreps);

		if ( !fork() ) //	child
		{
		String	command(getenv("DIESEL_EXE_DIR"));
			command += String("/refsel <refsel.in >refsel.out");
			execl("/bin/SH", "SH", "-c", command.chars(), NULL);
			cerr << "error executing " << command << endl;
			exit(1);
		}
		wait(0);


		iIrrep = irreps.first();
		while ( iIrrep )
		{
		char	buf[100];
		INT	irrep = irreps(iIrrep);
			cerr << "\tirrep=" << irrep << endl;
			sprintf(buf, "%d", irreps(iIrrep));
			chdir(buf);

		ifstream	is("refs.out");
		GrepAwk	ga(is);
		ConfigurationSet	refConfs, PTRefConfs;
			while ( !ga.illegal() )
			{
			String	buf(ga.getLine());
			istringstream	s(buf.chars(),istringstream::in | istringstream::out);
			Configuration<MOType>	conf(s);
				refConfs.add(conf);
				ga++;
			}

			cerr << "\t\treference space generation completed" << endl;
			writeSelectorInput("sel.in.all", 
				roots[irrep], nthExcitation, multiplicity, irrep, thresholds, 
				refConfs, PTRefConfs,
				activeReferenceThreshold!=0);


			chdir("..");
			irreps.next(iIrrep);
		}
	}




}


void	DieselInput::doDIESEL(INT multiplicity, INT irrep,
		const SLList<String> &thresholds,
		INT doNatOrb)
{

INT	selOK = 0;
	{
	ifstream	f("sel.out.all");
		if ( f )
		{
			selOK = 1;
		GrepAwk	ga(f);
			
		Pix	pix = thresholds.first();
			while ( pix )	
			{
				ga.head();
				ga.grep("R e s u l t s");
				ga.grep("threshold/mH");
				ga += 2;
			INT	ok = 0;
				while ( !ga.illegal() && ga.getNumberOfWords()>0 )
				{
					if ( fabs((atof(ga.getWord(1))/1000.0-
						atof(thresholds(pix)))/atof(thresholds(pix)))<1e-4 )
					{
						ok = 1;
						break;
					}
					ga++;
				}
				if ( !ok )
				{
					selOK = 0;
					break;
				}
				thresholds.next(pix);
			}			
		}
	}


	if ( selOK )
	{
		cerr << "\t\tusing previously performed selection" << endl;
	}
	else
	{
		cerr << "\t\tperforming selection on given thresholds" << endl;
	ifstream	selin("sel.in.all");
	GrepAwk	ga(selin);
	ofstream 	selin_o("sel.in.all");
		while ( !ga.illegal() )
		{
			if ( upcase(ga.getWord(1))=="SELECTIONTHRESHOLDS" )
			{
			Pix i;
				selin_o << "SelectionThresholds              = { "; 
				printSet(thresholds,selin_o); selin_o << "}" << endl;
			}
			else
				selin_o << ga.getLine() << endl;
			ga++;
		}
		if ( !fork() ) //	child
		{
		String	command(getenv("DIESEL_EXE_DIR"));
			command += String("/sel ") +
				(usePreviousSelection ? String("-r") : String("")) 
				+ String(" <sel.in.all >sel.out.all");
			execl("/bin/SH", "SH", "-c", command.chars(), NULL);
			cerr << "error executing " << command << endl;
			exit(1);
		}
		wait(0);
		if ( homogeneousSelection )
		{
			cerr << "\t\tcreating homogeneous selection" << endl;
		Pix	pix = thresholds.first();
			while ( pix )	
			{
				cerr << "\t\t\tthreshold " << thresholds(pix) << endl;
				if ( !fork() ) //	child
				{
				String	command(getenv("DIESEL_EXE_DIR"));
					command += String("/selsym ConfTree.dat.") + thresholds(pix)
						+ String(" -h >selsym.out.") + thresholds(pix);
					execl("/bin/SH", "SH", "-c", command.chars(), NULL);
					cerr << "error executing " << command << endl;
					exit(1);
				}
				wait(0);
				rename(String("ConfTree.dat.") + thresholds(pix) + String(".selsym"), 
					String("ConfTree.dat.") + thresholds(pix));
				
				thresholds.next(pix);
			}
			rename("sel.out.all", "sel.out.all.nohom");
			cerr << "\t\trecalculating energy sums" << endl;
			if ( !fork() ) //	child
			{
			String	command(getenv("DIESEL_EXE_DIR"));
				command += String("/sel -r") +
					+ String(" <sel.in.all >sel.out.all");
				execl("/bin/SH", "SH", "-c", command.chars(), NULL);
				cerr << "error executing " << command << endl;
				exit(1);
			}
			wait(0);
		
		}
	}

	writeDiagonalisatorInput("diag.in", RefThreshold, PTRefThreshold, roots[irrep]);

	cerr << "\t\tdiagonalization steps" << endl;
	{
	Pix	oldPix = NULL;
	Pix	pix = thresholds.first();
		while ( pix )	
		{
			cerr << "\t\t\tthreshold " << thresholds(pix);
	
	
//		ifstream	f(String("Eigenvectors.dat."+thresholds(pix)));
		ifstream	f(String("diag.out."+thresholds(pix)));
		INT	doit = 1;
			if ( f )
			{
			GrepAwk	ga(f);

				if ( ga.grep("all roots converged") )
					doit = 0;
			}
			if ( doit )
			{
				cerr << endl;
				unlink(String("ConfTree.dat"));
				symlink(String("ConfTree.dat.")+thresholds(pix), "ConfTree.dat");
				if ( !fork() ) //	child
				{
				String	command(getenv("DIESEL_EXE_DIR"));
					command += String("/diag -i <diag.in >diag.out.") + thresholds(pix);
					if ( oldPix )
					{
					FILE	*f = fopen(String("Eigenvectors.dat."+thresholds(oldPix)).chars(), "r");
						if ( f )
						{
							command += String(" -s ") + thresholds(oldPix);
							fclose(f);
						}
					}
					if ( nProc>1 )
					{
					char	b[1000];
						sprintf(b, " -p %d", nProc);
						command += b;
					}
					execl("/bin/SH", "SH", "-c", command.chars(), NULL);
					cerr << "error executing " << command << endl;
					exit(1);
				}
				wait(0);
				unlink(String("Eigenvectors.dat."+thresholds(pix)));
				rename("Eigenvectors.dat", String("Eigenvectors.dat."+thresholds(pix)));
			}
			else
			{
				cerr << ": already done" << endl;
			}
				
			oldPix = pix;
			thresholds.next(pix);
		}
	}
	cerr << "\t\tdiagonalization finished" << endl;



	if ( !doNatOrb )
	{
	INT	mrpt = 0;
	INT	mp3 = 0;
		{
		Pix	pix = fullMRCIExtrapolation.first();
			while ( pix )	
			{
				if ( upcase(fullMRCIExtrapolation(pix))=="MRMP2" )
					mrpt = 1;
				if ( upcase(fullMRCIExtrapolation(pix))=="MRMP3" )
				{
					mrpt = 1;
					mp3 = 1;
				}
				fullMRCIExtrapolation.next(pix);
			}
		}

		if ( mrpt )
		{
		ifstream	f(String("mrpt.out"));
		INT	doit = 1;
			if ( f )
			{
			GrepAwk	ga(f);

				if ( ga.grep("MR-MP2 Results:") )
					doit = 0;
				if ( mp3 )
				{
					doit = 1;
					if ( ga.grep("MR-MP2 Results:") )
						doit = 0;
				}
			}
			if ( doit )
			{
				cerr << "\t\tcalculating MRMP perturbation" << endl;
				writeMRPTInput("mrpt.in", roots[irrep], thresholds, mp3);

				if ( !fork() ) //	child
				{
				String	command(getenv("DIESEL_EXE_DIR"));
					command += String("/mrpt <mrpt.in >mrpt.out");
					execl("/bin/SH", "SH", "-c", command.chars(), NULL);
					cerr << "error executing " << command << endl;
					exit(1);
				}
				wait(0);
			}
			else
				cerr << "\t\tcalculating MRMP perturbation... already done" << endl;
			
		}
	}
}


void	DieselInput::createNatOrbs(Pix iMult, Pix iIrrep, INT &nNatOrbCalcs)
{
char	buf[100];

	fclose(fopen("=NatOrbRoot", "w"));
	cerr << "\t\tcreating natural orbitals" << endl;
	mkdir("NatOrb", S_IRUSR | S_IRGRP | S_IROTH | S_IXUSR | S_IXGRP | S_IXOTH | S_IWUSR);
	chdir("NatOrb");
	symlink(String(MOIntegralFile.name), "fort.31");
SLList<String> NatOrbThreshold;
		NatOrbThreshold.append(NaturalOrbitalSelectionThreshold);
INT irrep = irreps(iIrrep);
	doDIESEL(multiplicities(iMult), irrep, 
		NatOrbThreshold, 1);

	cerr << "\t\tcalculating density matrices" << endl;
	if ( !fork() ) //	child
	{
	String	command(getenv("DIESEL_EXE_DIR"));
		command += String("/dens -")+f31Format(MOIntegralFile.format)+String(" ")+MOLCASRootDir+String("motra.in ")
			+NaturalOrbitalSelectionThreshold+String(" all . all >dens.out");
		execl("/bin/SH", "SH", "-c", command.chars(), NULL);
		cerr << "error executing " << command << endl;
		exit(1);
	}
	wait(0);

	cerr << "\t\tcalculating natural orbitals" << endl;
	if ( averagedNaturalOrbitals )
	{
		nNatOrbCalcs = 1;
		if ( !fork() ) //	child
		{
		String	command(getenv("DIESEL_EXE_DIR"));
			command += "/natorb ";
			for ( unsigned INT i=0 ; i<roots[irrep].size() ; i++ )
			{
			char	buf[1000];
				sprintf(buf, "Density.dat.I.R%d_I.R%d.%s ",
					 i+1, i+1, NaturalOrbitalSelectionThreshold.chars());
				command += buf;
			}
			command += "-weight ";
			for ( unsigned INT i=0 ; i<roots[irrep].size() ; i++ )
				command += "1 ";

		char	buf[1000];
			sprintf(buf, "<%s%s >%sNATORB.%d",
				 MOLCASRootDir.chars(), orbitalFile.chars(), MOLCASRootDir.chars(), 1);
			command += buf;
			execl("/bin/SH", "SH", "-c", command.chars(), NULL);
			cerr << "error executing " << command << endl;
			exit(1);
		}
		wait(0);
	}
	else
	{
		nNatOrbCalcs = roots.size();
		for ( INT i=0 ; i<nNatOrbCalcs ; i++ )
		{
			if ( !fork() ) //	child
			{
			String	command(getenv("DIESEL_EXE_DIR"));
			char	buf[1000];
				sprintf(buf, "Density.dat.I.R%d_I.R%d.%s <%s%s >%sNATORB.%d",
					 i+1, i+1, NaturalOrbitalSelectionThreshold.chars(),
					 MOLCASRootDir.chars(), orbitalFile.chars(), MOLCASRootDir.chars(), i+1);
				command += String("/natorb ")+String(buf);
				execl("/bin/SH", "SH", "-c", command.chars(), NULL);
				cerr << "error executing " << command << endl;
				exit(1);
			}
			wait(0);
		}
	}
char	current[1000];

	chdir("..");
	getcwd(current, 1000);
	chdir(MOLCASRootDir.chars());

	rename("STONEY", "STONEY.BASE");
	for ( INT i=0 ; i<nNatOrbCalcs ; i++ )
	{
		cerr << "\t\tperforming MO transformation for root #" << i+1 << endl;
		unlink(String("INPORB"));
	char	buf[10000];
		sprintf(buf, "%d", i+1);
		symlink(String("NATORB.")+String(buf), "INPORB");
		if ( !fork() ) //	child
		{
		String	command(getenv("MOLCAS_EXE_DIR"));
			command += String("/motra <motra.in >motra.out");
			execl("/bin/SH", "SH", "-c", command.chars(), NULL);
			cerr << "error executing " << command << endl;
			exit(1);
		}
		wait(0);
		if ( !fork() ) //	child
		{
		String	command(getenv("MOLCAS_EXE_DIR"));
			command += String("/form31.small <form31.in >form31.out");
			execl("/bin/SH", "SH", "-c", command.chars(), NULL);
			cerr << "error executing " << command << endl;
			exit(1);
		}
		wait(0);
		unlink(String("STONEY.NatOrb.")+String(buf));
		rename("STONEY", String("STONEY.NatOrb.")+String(buf));
	}
	rename("STONEY.BASE", "STONEY");


	chdir(current);
	for ( INT i=0 ; i<nNatOrbCalcs ; i++ )
	{
		cerr << "\t\tperforming MR-CI calculation with natural orbitals for root #" << i+1 << endl;
		sprintf(buf, "%d", i+1);
		mkdir(buf, S_IRUSR | S_IRGRP | S_IROTH | S_IXUSR | S_IXGRP | S_IXOTH | S_IWUSR);
		chdir(buf);

		sprintf(buf, "%sSTONEY.NatOrb.%d", MOLCASRootDir.chars(), i+1);
		symlink(String(buf), "fort.31");
		doDIESEL(multiplicities(iMult), irreps(iIrrep), thresholds, 0);
		chdir("..");
	}
}



void	DieselInput::start()
{
	cerr << "*******************************************************************************" << endl;
	cerr << "*                                                                             *" << endl;
	cerr << "*                                diesel protocol                              *" << endl;
	cerr << "*                                                                             *" << endl;
	cerr << "*******************************************************************************" << endl;


	cerr << endl;
	cerr << endl;
	fclose(fopen("=Multiplicity", "w"));
	
	mkdir(TempDir.chars(), S_IRUSR | S_IRGRP | S_IROTH | S_IXUSR | S_IXGRP | S_IXOTH | S_IWUSR);

char	buf[100];
Pix	iMult = multiplicities.first();
	while ( iMult )
	{
		cerr << "multiplicity=" << multiplicities(iMult) << endl;
		sprintf(buf, "%d", multiplicities(iMult));
		mkdir(buf, S_IRUSR | S_IRGRP | S_IROTH | S_IXUSR | S_IXGRP | S_IXOTH | S_IWUSR);
		chdir(buf);
		fclose(fopen("=Irrep", "w"));
//		symlink(String(MOIntegralFile.name), "fort.31");
		createReferenceSpace(multiplicities(iMult));



	Pix	iIrrep = irreps.first();
		while ( iIrrep )
		{
		INT irrep = irreps(iIrrep);
			cerr << "\tirrep=" << irrep << endl;
			sprintf(buf, "%d", irreps(iIrrep));
			mkdir(buf, S_IRUSR | S_IRGRP | S_IROTH | S_IXUSR | S_IXGRP | S_IXOTH | S_IWUSR);
			chdir(buf);
			
		INT	nNatOrbCalcs = 1;
			if ( useNaturalOrbitals )
				createNatOrbs(iMult, iIrrep, nNatOrbCalcs);
			else
			{
				symlink(String(MOIntegralFile.name), "fort.31");
				doDIESEL(multiplicities(iMult), irreps(iIrrep), thresholds, 0);
			}
		
			if ( !thresholds.length() )
				goto multcont;
		
			for ( INT i=0 ; i<nNatOrbCalcs ; i++ )
			{
//			String	outName("../../diesel.out");
				if ( useNaturalOrbitals )
				{
				char	buf[100];
					sprintf(buf, "%d", i+1);
					chdir(buf);
//					outName = "../" + outName;
				}
				cout << "============================================" << endl;
				{
				char	command[1000];
					sprintf(command, "%s/getDirInfo", getenv("DIESEL_EXE_DIR"));
					system(command);
				}
				{
					cout << "============================================" << endl;
					cout << endl;
					cout << endl;
					cout << endl;


				ifstream	seloutall("sel.out.all");
				GrepAwk	ga(seloutall);

					if ( ga.grep("selected reference configuration") )
					{
						cout << ga.getLine() << endl;
						ga++;
					INT ii = 0;
						while ( !ga.illegal() && ga.getNumberOfWords()>0 )
						{
							cout << ga.getLine() << "      # " << ++ii << endl;
							ga++;
						}
					}
					else
					{
						ga.head();
						ga.grep("RefConfs");
						ga += 2;

						cout << "reference configurations:" << endl;
						while ( !ga.illegal() && ga.getWord(1)!="}" )
						{
							cout << ga.getLine() << endl;
							ga++;
						}
					}
					cout << endl;
					cout << endl;

					ga.head();
					if ( !ga.grep("selected reference configuration") )
						ga.head();
					ga.grep("number of configurations:");
					cout << ga.getLine() << endl;
					ga++;
					cout << ga.getLine() << endl;
					ga++;
					cout << endl;
					cout << endl;
					cout << "selected roots within reference space:" << endl;
					ga.grep("eigenvalues of reference matrix:");
				INT	r = 0;
					while ( !ga.illegal() && ga.getNumberOfWords()>0 )
					{
						if ( ga.getWord(2)=="<---" )
							cout << "root #" << ++r << ": " << ga.getWord(1) << endl;
						ga++;
					}

					if ( verbosity.isActive(Verbosity::RefMatEigenVectors) )
					{
						cout.setf(ios::fixed, ios::floatfield);
						cout << endl;
						cout << endl;
						cout << "eigenvectors of reference matrix (columnwise):" << endl;
						ga.grep("eigenvectors of reference matrix (columnwise):");
						ga++;
						ga++;
					INT	 dim = 3*getRoots(irrep).size()/2;
						cout << setw(3) << "ref" << " ";
						for ( INT j=0 ; j<dim ; j++ )
							cout << setw(7) << setprecision(0) << j+1 << " ";
						cout << endl;
						while ( !ga.illegal() && ga.getNumberOfWords()>0 )
						{
							cout << setw(3) << ga.getWord(1) << " ";
							for ( INT i=2 ; i<=1+dim ; i++ )
							{
							double	h = atof(ga.getWord(i));
								cout << setw(7) << setprecision(4) << h << " ";
							}
							cout << endl;
							ga++;
						}


						cout << endl;
						cout << endl;
					}


					cout << endl;
					cout << endl;
					ga.grep("total generated configurations/CSFs:");
					for ( INT i=0 ; i<8 ; i++ )
					{
						cout << ga.getLine() << endl;
						ga++;
					}
					cout << endl;
					cout << endl;
					cout << endl;

					{
					char	command[1000];
						sprintf(command, "%s/dr", getenv("DIESEL_EXE_DIR"));
						system(command);
					}
				char newpage = 12;
					cout << newpage;



				}
				if ( useNaturalOrbitals )
					chdir("..");
			}
			


multcont:
			chdir(".."),
			irreps.next(iIrrep);
			cerr << endl;
		}

		if ( propertyThresholds.length() )
		{
			cerr << "\tcalculating one particle density matrices" << endl;


		Pix	iPropThresh = propertyThresholds.first();
			while ( iPropThresh )
			{
				cerr << "\t\tthreshold=" << propertyThresholds(iPropThresh) << endl;
				symlink(String(MOIntegralFile.name), "fort.31");

			Pix	iIrrep = irreps.first();
				while ( iIrrep )
				{
					cerr << "\t\t\tiirrep=" << irreps(iIrrep) << flush;

				char	buf[1000];
					sprintf(buf, "Densities/Density.dat.I%dR%d_I%dR%d.%s",
						 irreps(iIrrep), 1,
						 irreps(iIrrep), 1,
						 propertyThresholds(iPropThresh).chars());
					{
					ifstream	fdens(buf);
						if ( !fdens )
						{
						char	command[1000];
							sprintf(command, "%s/dens -%c %smotra.in %s all %d all >>dens.out.%s",
								getenv("DIESEL_EXE_DIR"), 
								f31Format(MOIntegralFile.format), MOLCASRootDir.chars(),
								propertyThresholds(iPropThresh).chars(),
								irreps(iIrrep),
								propertyThresholds(iPropThresh).chars());
							system(command);
						}
						else
							cerr << "  already done";
					}


					cerr << endl;

				Pix	jIrrep = iIrrep;
					irreps.next(jIrrep);
					while ( jIrrep )
					{
						cerr << "\t\t\t\tjirrep=" << irreps(jIrrep) << flush;

				char	buf[1000];
					sprintf(buf, "Densities/Density.dat.I%dR%d_I%dR%d.%s",
						 irreps(iIrrep), 1,
						 irreps(jIrrep), 1,
						 propertyThresholds(iPropThresh).chars());
						{
						ifstream	fdens(buf);
							if ( !fdens )
							{
						char	command[1000];
							sprintf(command, "%s/dens -%c %smotra.in %s all %d all %d >>dens.out.%s",
								getenv("DIESEL_EXE_DIR"), 
								f31Format(MOIntegralFile.format), MOLCASRootDir.chars(),
								propertyThresholds(iPropThresh).chars(),
								irreps(iIrrep), irreps(jIrrep),
								propertyThresholds(iPropThresh).chars());
								system(command);
							}
							else
								cerr << "  already done";
						}
						cerr << endl;

						irreps.next(jIrrep);
					}

					irreps.next(iIrrep);
				}



			ifstream	propOutI("prop.dat."+propertyThresholds(iPropThresh));
				if ( !propOutI || true )
				{
// Bug fix for NFS mounted filesystems:
// do no longer rename Density*.dat files due to probable cache incoherencies
// 03.11.2000
/*					{
					DIR	*dir = opendir(".");
					struct dirent	*entry;

						while ( (entry=readdir(dir)) )
						{
							if ( !strncmp("Density.dat.",entry->d_name, 12)  )
								rename(entry->d_name, 
									String("Densities/"+String(entry->d_name)).chars());
						}
						closedir(dir);
					}
*/
					cerr << "\t\t\tcalculating properties" << endl;

				SLList<String>	ls;
					{
						DIR	*dir = opendir("Densities");
						struct dirent	*entry;

							while ( (entry=readdir(dir)) )
							{
								if ( !strncmp("Density.dat.", entry->d_name, 12) &&
									propertyThresholds(iPropThresh)==
										entry->d_name+strlen(entry->d_name)-propertyThresholds(iPropThresh).length() )

								{
								char	d[1000];
									sprintf(d, "Densities/%s", entry->d_name);
								String	S(d);
									ls.append(S);
								}
							}
							closedir(dir);
					}


					//	cout << args << endl;


				pid_t	pid = fork();
					if ( !pid ) //	child
					{
					char	buf3[1000];
						sprintf(buf3, "prop.dat.%s", propertyThresholds(iPropThresh).chars());
					INT	f = open(buf3, O_CREAT | O_TRUNC | O_WRONLY, S_IRUSR | S_IWUSR);
						close(1);
						dup(f);
					String	command(getenv("DIESEL_EXE_DIR"));
						command += "/prop";
					const char	*argsV[ls.length()+5+1];
						argsV[0] = "prop";
						argsV[1] = "MLTPL1";
					char	buf1[1000];
						sprintf(buf1, "%s%s",
									MOLCASRootDir.chars(), orbitalFile.chars());
						argsV[2] = buf1;
					char	buf2[1000];
						sprintf(buf2, "%sONEINT", MOLCASRootDir.chars());
						argsV[3] = buf2;
//					char	buf3[1000];
//						sprintf(buf3, ">prop.dat.%s", propertyThresholds(iPropThresh).chars());
//						argsV[3] = buf3;
					Pix pix = ls.first();
					INT	i = 4;
						while ( pix )
						{
							argsV[i++] = ls(pix);
							ls.next(pix);
						}
						argsV[ls.length()+4] = NULL;
						for ( INT ii=0 ; ; ++ii )
						{
							if ( !argsV[ii] )
								break;
						}

						execv(command.chars(), (char *const *) argsV);

					}
					wait(0);
					cerr << "\t\t\tproperty calculation finished" << endl;
				}
				else
					cerr << "\t\t\tproperty calculation already done" << endl;

				propertyThresholds.next(iPropThresh);
			}
				cerr << "\tone particle density matrices calculation finished" << endl;

		}

		chdir(".."),
		multiplicities.next(iMult);
		cerr << endl;
	}
	
	rmdir(TempDir.chars());
	
}









void	DieselInput::setMORestriction(const char *yytext)
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
	moRestriction = new MORestriction(s);
}


void	DieselInput::setMOEquivalence(const char *yytext, INT auto1)
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

#include "../../../../lib/Container/SLList.cc"

template class SLList<String>;
template class SLList<double>;

