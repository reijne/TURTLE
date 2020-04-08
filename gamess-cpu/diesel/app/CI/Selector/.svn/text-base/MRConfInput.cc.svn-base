//***********************************************************************
//
//	Name:			MRConfInput.cc
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

#include "MRConfInput.h"

#include <stdio.h>
#include <stdlib.h>
#include <string>

#include <sstream>

#include "../../../lib/QM/MO/MRMOs.h"
#include "../../../lib/QM/MO/MOMapping.h"

#include "../../../lib/QM/IntegralContainer/DiagIntContainer.h"

#include "../../../lib/QM/IO/Verbosity.h"

#include "../../../lib/QM/MO/Iterators/MORestriction.h"
#include "../../../lib/QM/MO/Iterators/MOEquivalence.h"
#include "../../../lib/QM/MO/MOStatistics.h"

#include "../../../lib/Container/AVLSet.h"
#include "../../../lib/QM/MO/Iterators/MOIterator.h"

using namespace std;


MRConfInput::MRConfInput()
{
	NumberOfElectrons = 0;
	Multiplicity = 0;
	maxRefOpenShells = 0;
	firstGuessConfs = 0;
	selectInternal = 0;
	selectNthExcitation = 0;
	storePTEnergy = 0;
	storePTCoef = 0;
	activeSpaceExcitationLevel = 0;
	autoRef = 0;
	autoEquiv = 0;
	irrep = 0;
	estimationMode = EnergyMap::EpsteinNesbet;
	nSelectionThresholds = 0;
	SelectionThreshold = NULL;
	SelectionThresholdString = NULL;
	randomProb = 0;
	NumberOfRoots = 0;
	rootNumbers = NULL;
	rootEnergies = NULL;
	moRestriction = NULL;
	moEquivalence = NULL;
	moStatistics = NULL;
	mrmos = 0;
	activeReferenceThreshold = 0;
}


MRConfInput::~MRConfInput()
{
	if ( mrmos )
		delete mrmos;
		
	if ( rootNumbers )
		free(rootNumbers);

	if ( rootEnergies )
		free(rootEnergies);
		
	if ( SelectionThreshold )
		free(SelectionThreshold);

	if ( SelectionThresholdString )
	{
		for ( INT i=0 ; i<nSelectionThresholds ; i++ )
			delete[] SelectionThresholdString[i];
		free(SelectionThresholdString);
	}

	if ( moRestriction )
		delete moRestriction;

	if ( moEquivalence )
		delete moEquivalence;

	if ( moStatistics )
		delete moStatistics;
}

MRConfInput::MRConfInput(const MRConfInput &mrconf) :
	Refs(((MRConfInput *) &mrconf)->getRefConfSet()), 
	PTRefs(((MRConfInput *) &mrconf)->getPTRefConfSet())
{
	MOIntegralFile = mrconf.getMOIntegralFile();
	NumberOfElectrons = mrconf.getNumberOfElectrons();
	ExcitationLevel = mrconf.getExcitationLevel();
	Multiplicity = mrconf.getMultiplicity();
	maxRefOpenShells = mrconf.maxRefOpenShells;
	firstGuessConfs = mrconf.firstGuessConfs;
	selectInternal = mrconf.getSelectInternal();
	selectNthExcitation = mrconf.getSelectNthExcitation();
	storePTEnergy = mrconf.getStorePTEnergy();
	storePTCoef = mrconf.getStorePTCoef();
	irrep = mrconf.getIrRep();
	autoRef = mrconf.isAutoRef();
	autoEquiv = mrconf.isAutoEquiv();
	activeSpaceExcitationLevel = mrconf.activeSpaceExcitationLevel;
	activeReferenceThreshold = mrconf.activeReferenceThreshold;
	randomProb = mrconf.randomProb;

	estimationMode = mrconf.estimationMode;

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

	annihilatorSpace = mrconf.annihilatorSpace;
	creatorSpace = mrconf.creatorSpace;

	NumberOfRoots = mrconf.getNumberOfRoots();
	rootNumbers = (INT *) malloc(NumberOfRoots * sizeof(INT));
	for ( INT i=0 ; i<NumberOfRoots ; i++ )
		rootNumbers[i] = mrconf.getRootNumber(i);
	mrmos = new MRMOs(*mrconf.getMRMOs());
	
	if ( mrconf.getMORestriction() )
		moRestriction = new MORestriction(*mrconf.getMORestriction());
	else
		moRestriction = NULL;

	if ( mrconf.getMOEquivalence() )
		moEquivalence = new MOEquivalence(*mrconf.getMOEquivalence());
	else
		moEquivalence = NULL;

	if ( mrconf.getMOStatistics() )
		moStatistics = new MOStatistics(*mrconf.getMOStatistics());
	else
		moStatistics = NULL;


	if ( mrconf.rootEnergies )
	{
		rootEnergies = (EnergyType *) malloc(NumberOfRoots * sizeof(EnergyType));
		for ( INT i=0 ; i<NumberOfRoots ; i++ )
			rootEnergies[i] = mrconf.getRootEnergy(i);
	}
	else
		rootEnergies = NULL;
		
	moMappingNoPTRef = mrconf.moMappingNoPTRef;
}

MRConfInput & MRConfInput::operator = (const MRConfInput &mrconf)
{
	Refs = ConfigurationSet(((MRConfInput *) &mrconf)->getRefConfSet());
	PTRefs = ConfigurationSet(((MRConfInput *) &mrconf)->getPTRefConfSet());
	
	MOIntegralFile = mrconf.getMOIntegralFile();
	NumberOfElectrons = mrconf.getNumberOfElectrons();
	ExcitationLevel = mrconf.getExcitationLevel();
	Multiplicity = mrconf.getMultiplicity();
	maxRefOpenShells = mrconf.maxRefOpenShells;
	firstGuessConfs = mrconf.firstGuessConfs;
	selectInternal = mrconf.getSelectInternal();
	selectNthExcitation = mrconf.getSelectNthExcitation();
	storePTEnergy = mrconf.getStorePTEnergy();
	storePTCoef = mrconf.getStorePTCoef();
	irrep = mrconf.getIrRep();
	autoRef = mrconf.isAutoRef();
	autoEquiv = mrconf.isAutoEquiv();
	activeSpaceExcitationLevel = mrconf.activeSpaceExcitationLevel;
	activeReferenceThreshold = mrconf.activeReferenceThreshold;

	randomProb = mrconf.randomProb;

	estimationMode = mrconf.estimationMode;

	annihilatorSpace = mrconf.annihilatorSpace;
	creatorSpace = mrconf.creatorSpace;

	if ( mrmos )
		delete mrmos;
		
	if ( rootNumbers )
		free(rootNumbers);

	if ( rootEnergies )
		free(rootEnergies);
		
	if ( SelectionThreshold )
		free(SelectionThreshold);

	if ( SelectionThresholdString )
	{
		for ( INT i=0 ; i<nSelectionThresholds ; i++ )
			delete SelectionThresholdString[i];
		free(SelectionThresholdString);
	}

	if ( moRestriction )
		delete moRestriction;

	if ( moEquivalence )
		delete moEquivalence;

	if ( moStatistics )
		delete moStatistics;


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
	
	if ( mrconf.getMORestriction() )
		moRestriction = new MORestriction(*mrconf.getMORestriction());
	else
		moRestriction = NULL;

	if ( mrconf.getMOEquivalence() )
		moEquivalence = new MOEquivalence(*mrconf.getMOEquivalence());
	else
		moEquivalence = NULL;

	if ( mrconf.getMOStatistics() )
		moStatistics = new MOStatistics(*mrconf.getMOStatistics());
	else
		moStatistics = NULL;


	if ( mrconf.rootEnergies )
	{
		rootEnergies = (EnergyType *) malloc(NumberOfRoots * sizeof(EnergyType));
		for ( INT i=0 ; i<NumberOfRoots ; i++ )
			rootEnergies[i] = mrconf.getRootEnergy(i);
	}
	else
		rootEnergies = NULL;

	moMappingNoPTRef = mrconf.moMappingNoPTRef;

	return *this;
}

#define printSet(x) i = x.first(); while ( i ) { s << (x)(i) << " "; (x).next(i); }

ostream& operator<<(ostream& s, const MRConfInput & mrconf)
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
	s << "EstimationMode             = ";
	switch ( mrconf.estimationMode )
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
	if ( mrconf.moRestriction )
		s << "MORestrictions             = {" << mrconf.moRestriction->getRestriction() << "}" << endl;
	else
		s << "MORestrictions             = none" << endl;

	if ( mrconf.moEquivalence )
		s << "MOEquivalence              = {" << mrconf.moEquivalence->getEquivalence() << "}" << endl;
	else
		s << "MOEquivalence              = none" << endl;

	s << "MOStatistics               = " << (mrconf.moStatistics ? "yes" : "no") << endl;

	s << "NumberOfElectrons          = " << mrconf.NumberOfElectrons << endl;
	s << "ExcitationLevel            = " << mrconf.ExcitationLevel << endl;
	s << "Multiplicity               = " << mrconf.Multiplicity << endl;
	s << "IrRep                      = " << mrconf.irrep << endl;

	s << "Roots                      = { ";
	for ( INT i=0 ; i<mrconf.NumberOfRoots ; i++ )
		s << mrconf.getRootNumber(i) << " ";
	s << "}" << endl;

	if ( mrconf.rootEnergies )
	{
		s << "RootEnergies               = { ";
		for ( INT i=0 ; i<mrconf.NumberOfRoots ; i++ )
			s << mrconf.getRootEnergy(i) << " ";
		s << "}" << endl;
	}
	else
		s << "RootEnergies               = reference eigenvalues" << endl;

	s << "SelectionThresholds        = { ";
	for ( INT i=0 ; i<mrconf.getNumberOfSelectionThresholds() ; i++ )
		s << mrconf.getSelectionThreshold(i) << " ";
	s << "}" << endl;

	s << "RandomProb                 = ";	s << 100.0*mrconf.getRandomProb() << " ";
	s << " %" << endl;


	{
	Pix	i = NULL;
		if ( mrconf.annihilatorSpace.length() )
		{
			s << "AnnihilatorSpace           = { "; 
			printSet(mrconf.annihilatorSpace); 
			s << "}" << endl;
		}
		else
			s << "AnnihilatorSpace           = inactive" << endl;

		if ( mrconf.creatorSpace.length() )
		{
			s << "CreatorSpace               = { "; 
			printSet(mrconf.creatorSpace); 
			s << "}" << endl;
		}
		else
			s << "CreatorSpace               = inactive" << endl;
	}

	s << "ActiveSpaceExcitationLevel = " << mrconf.activeSpaceExcitationLevel << endl;
	s << "maxRefOpenShells           = " << mrconf.maxRefOpenShells << endl;
	s << "activeReferenceThreshold   = " << mrconf.activeReferenceThreshold << endl;


	s << "FirstGuessConfs            = " << mrconf.firstGuessConfs << endl;
	s << "selectInternal             = " << ( mrconf.selectInternal ? "yes" : "no") << endl;
	s << "selectNthExcitation        = { ";
	for ( INT i=0 ; i<(INT) sizeof(INT)*8 ; i++ )
		if ( mrconf.selectNthExcitation & (1<<i) )
			s << i+1 << " ";
	s << "}" << endl;

	s << "StorePTEnergy              = " << ( mrconf.storePTEnergy ? "yes" : "no") << endl;
	s << "StorePTCoef                = " << ( mrconf.storePTCoef ? "yes" : "no") << endl;

	s << "RefConfs                   = " // << mrconf.Refs.n
		<< endl;
	s << "{" << endl;
	((MRConfInput *) &mrconf)->Refs.setNumbering(1);
	s << mrconf.Refs;
	((MRConfInput *) &mrconf)->Refs.setNumbering(0);
	s << "}" << endl;

	s << "PTRefConfs                 = ";
	if ( mrconf.PTRefs==mrconf.Refs )
		cout << "RefConfs" << endl;
	else
	{
		cout << endl;
		s << "{" << endl;
		((MRConfInput *) &mrconf)->PTRefs.setNumbering(1);
		s << mrconf.PTRefs;
		((MRConfInput *) &mrconf)->PTRefs.setNumbering(0);
		s << "}" << endl;
	}
	return s;
}



void	MRConfInput::initInternalExternal()
{
INT internal[mrmos->getMaxMO()];
	memset(internal, 0, mrmos->getMaxMO()*sizeof(INT));

	{
	Pix	i = Refs.first();
	
		while ( i )
		{
			for ( INT j=0 ; j<getRefConfSet()(i).getNumberOfOpenShells() ; j++ )
				internal[getRefConfSet()(i).getOpenShell(j) - 1] = 1;
			for ( INT j=0 ; j<getRefConfSet()(i).getNumberOfClosedShells() ; j++ )
				internal[getRefConfSet()(i).getClosedShell(j) - 1] = 1;
				
			Refs.next(i);
		}
	}
	getMOMappingNoPTRef() = MOMapping(mrmos->getMaxMO(), internal, *mrmos);
	
	{
	Pix	i = PTRefs.first();
	
		while ( i )
		{
			for ( INT j=0 ; j<getPTRefConfSet()(i).getNumberOfOpenShells() ; j++ )
				internal[getPTRefConfSet()(i).getOpenShell(j) - 1] = 1;
			for ( INT j=0 ; j<getPTRefConfSet()(i).getNumberOfClosedShells() ; j++ )
				internal[getPTRefConfSet()(i).getClosedShell(j) - 1] = 1;
				
			PTRefs.next(i);
		}
	}
	
	moMapping = MOMapping(mrmos->getMaxMO(), internal, *mrmos);


	// apply MO-mapping to reference configurations
Configuration<MOType>	inactive;
	{
	Pix	i = Refs.first();
		if ( i )
			inactive = getRefConfSet()(i);
	
		while ( i )
		{
			for ( INT j=0 ; j<getRefConfSet()(i).getNumberOfOpenShells() ; j++ )
			{
				getRefConfSet()(i).setOpenShell(j, moMapping.getContinuous(
					getRefConfSet()(i).getOpenShell(j)));
				mrmos->setInternal(getRefConfSet()(i).getOpenShell(j));
			}
			for ( INT j=0 ; j<getRefConfSet()(i).getNumberOfClosedShells() ; j++ )
			{
				getRefConfSet()(i).setClosedShell(j, moMapping.getContinuous(
					getRefConfSet()(i).getClosedShell(j)));
				mrmos->setInternal(getRefConfSet()(i).getClosedShell(j));
			}
			inactive &= getRefConfSet()(i);
			Refs.next(i);
		}
	}

	// apply MO-mapping to PT reference configurations
	{
	Pix	i = PTRefs.first();
	
		while ( i )
		{
			for ( INT j=0 ; j<getPTRefConfSet()(i).getNumberOfOpenShells() ; j++ )
			{
				getPTRefConfSet()(i).setOpenShell(j, moMapping.getContinuous(
					getPTRefConfSet()(i).getOpenShell(j)));
				mrmos->setInternal(getPTRefConfSet()(i).getOpenShell(j));
			}
			for ( INT j=0 ; j<getPTRefConfSet()(i).getNumberOfClosedShells() ; j++ )
			{
				getPTRefConfSet()(i).setClosedShell(j, moMapping.getContinuous(
					getPTRefConfSet()(i).getClosedShell(j)));
				mrmos->setInternal(getPTRefConfSet()(i).getClosedShell(j));
			}
			PTRefs.next(i);
		}
	}
	mrmos->setNumberOfInactiveMOs(
		0*inactive.getNumberOfOpenShells() + inactive.getNumberOfClosedShells());

/*	for ( INT i=0 ; i<inactive.getNumberOfOpenShells() ; i++ )
		mrmos->setInactiveMO(i, inactive.getOpenShell(i));
*/
	for ( INT i=0 ; i<inactive.getNumberOfClosedShells() ; i++ )
		mrmos->setInactiveMO(0*inactive.getNumberOfOpenShells() + i,
			 inactive.getClosedShell(i));

	if ( verbosity.isActive(Verbosity::MOs) )
		cout << moMapping << endl;

	mrmos->initIntExt();

	if ( verbosity.isActive(Verbosity::MOs) )
	{
		cout << *mrmos << endl;
		cout << "no errors." << endl << endl;
	}
}


#include <FlexLexer.h>

MRConfInput *LexMRConf;
INT	LexMRConfErrors;



static void	ConfError(INT i, Configuration<MOType> conf)
{
	cout << "Error(s) at conf " << i << ": \"" << conf << "\":" << endl;
}


INT	MRConfInput::checkConfs(ConfigurationSet &confSet,
	INT &els, INT &elsErrors, INT &errors)
{
Pix	i = confSet.first();
INT	ii = 1;
	while ( i )
	{
ostringstream	errString(ostringstream::in | ostringstream::out);
		errString.clear();

                streampos firstpos = errString.tellp();

		if ( els == -1 )
			els = confSet(i).getNumberOfElectrons();
		else
		{
			if ( els != confSet(i).getNumberOfElectrons() )
			{
				elsErrors++;
				if ( NumberOfElectrons == -1 )
					errString << "inconsistent number of electrons: "
						<< confSet(i).getNumberOfElectrons() << endl;
				else
					errString << "wrong number of electrons: "
						<< confSet(i).getNumberOfElectrons() 
						<< " should be " << els << endl;
			}
			else
				if ( elsErrors && NumberOfElectrons == -1 )
					errString << "inconsistent number of electrons: "
						<< confSet(i).getNumberOfElectrons() << endl;
		}

		if ( Multiplicity>0 )
		{
			if ( (confSet(i).getNumberOfOpenShells() & 1) ==
					(Multiplicity & 1) )
			{
				errString << "number of open shells must be " 
					<<  ((Multiplicity & 1) ? "even" : "odd") << "." << endl;
				errors++;
			}

			if ( confSet(i).getNumberOfOpenShells()<Multiplicity-1 )
			{
				errString << "number of open shells too small." << endl;
				errors++;
			}
		}

		if ( confSet(i).check() )
		{
			errString << "MO numbering is wrong." << endl;
			errors++;
		}
				
		for ( INT j=0 ; j<confSet(i).getNumberOfOpenShells() ; j++ )
		{
			if ( confSet(i).getOpenShell(j)>mrmos->getMaxMO() )
			{
				errString << "MO number " 
					<< confSet(i).getOpenShell(j) << " too large." << endl;
				errors++;
			}
			if ( confSet(i).getOpenShell(j)<=0 )
			{
				errString << "MO number " 
					<< confSet(i).getOpenShell(j) << " too small." << endl;
				errors++;
			}
		}
		for ( INT j=0 ; j<confSet(i).getNumberOfClosedShells() ; j++ )
		{
			if ( confSet(i).getClosedShell(j)>mrmos->getMaxMO() )
			{
				errString << "MO number " 
					<< confSet(i).getClosedShell(j) << " too large." << endl;
				errors++;
			}
			if ( confSet(i).getClosedShell(j)<=0 )
			{
				errString << "MO number " 
					<< confSet(i).getClosedShell(j) << " too small." << endl;
				errors++;
			}
		}
		
		{
		IrRep	Sym = 0;
			for ( INT j=0 ; j<confSet(i).getNumberOfOpenShells() ; j++ )
				Sym = mrmos->getProd(Sym, 
					mrmos->getIrRep(confSet(i).getOpenShell(j)));

			if ( irrep == -1 )
				irrep = Sym;
				
			if ( Sym != irrep )
			{
				errString << "symmetry is wrong: "
					<< Sym << " should be " << irrep << "." << endl;
				errors++;
			}
		}
			
		

		if ( errString.tellp() != firstpos)
		{
			ConfError(ii, confSet(i));
			cout << errString.str() << endl;
		}
		
		
		confSet.next(i);
		ii++;
	}
	return errors;
}



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

istream& operator>>(istream& s, MRConfInput & mrconf)
{
yyFlexLexer	lexer(&s);

ostream	*outStream = new ostringstream(ostringstream::in | ostringstream::out);
	LexMRConfErrors = 0;

	mrconf.NumberOfElectrons = -1;
	mrconf.Refs.clear();
	mrconf.PTRefs.clear();
	mrconf.NumberOfRoots = 0;
	mrconf.NumberOfRootsE = 0;
	mrconf.rootNumbers = NULL;
	mrconf.rootEnergies = NULL;
	mrconf.selectInternal = -1;
	mrconf.ExcitationLevel = -1;
	mrconf.storePTEnergy = -1;
	mrconf.storePTCoef = -1;
	mrconf.Multiplicity = -1;
	mrconf.maxRefOpenShells = -1;
	mrconf.firstGuessConfs = -1;
	mrconf.activeSpaceExcitationLevel = -1;
	mrconf.activeReferenceThreshold = -1;
	mrconf.randomProb = -1;
	mrconf.autoRef = 0;
	mrconf.autoEquiv = 0;
	mrconf.irrep = -1;
	mrconf.estimationMode = EnergyMap::EpsteinNesbet;
	mrconf.nSelectionThresholds = 0;
	mrconf.SelectionThreshold = NULL;
	mrconf.SelectionThresholdString = NULL;
	mrconf.moRestriction = NULL;
	mrconf.moEquivalence = NULL;
	mrconf.moStatistics = NULL;
	mrconf.MOIntegralFile = Fort31File("", Fort31RecordFormatAuto);
	
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

	if ( mrconf.NumberOfElectrons == -1 )
	{
		if ( mrconf.autoRef )
		{
			cout << "\"NumberOfElectrons\" is mandatory when using \"RefConfs = auto\"" << endl;
			errors++;
		}
		else
			*outStream << "\"NumberOfElectrons\" is missing, defaults to configuration input." << endl;
	}

	if ( mrconf.activeSpaceExcitationLevel == -1 )
	{
		*outStream << "\"ActiveSpaceExcitationLevel\" is missing, defaults to 1." << endl;
		mrconf.activeSpaceExcitationLevel = 1;
	}

	if ( mrconf.randomProb == -1 )
	{
		*outStream << "\"RandomProb\" is missing, defaults to 0 % (inactive)." << endl;
		mrconf.randomProb = 0;
	}

	if ( mrconf.ExcitationLevel == -1 )
	{
		*outStream << "\"ExcitationLevel\" is missing, defaults to 2." << endl;
		mrconf.ExcitationLevel = 2;
	}

	if ( mrconf.selectInternal == -1 )
	{
		*outStream << "\"SelectInternal\" is missing, defaults to \"no\"." << endl;
		mrconf.selectInternal = 0;
	}

	if ( mrconf.storePTEnergy == -1 )
	{
		*outStream << "\"StorePTEnergy\" is missing, defaults to \"no\"." << endl;
		mrconf.storePTEnergy = 0;
	}

	if ( mrconf.storePTCoef == -1 )
	{
		*outStream << "\"StorePTCoef\" is missing, defaults to \"no\"." << endl;
		mrconf.storePTCoef = 0;
	}

	if ( mrconf.irrep == -1 )
	{
		if ( mrconf.autoRef )
		{
			cout << "\"IrRep\" is mandatory when using \"RefConfs = auto\"" << endl;
			errors++;
		}
		*outStream << "\"IrRep\" is missing, defaults to configuration input." << endl;	
	}

	if ( strlen(mrconf.MOIntegralFile.name) == 0 )
	{
		*outStream << "\"MOIntegralFilename\" is missing, defaults to \"fort.31\"." << endl;	
		mrconf.MOIntegralFile.name = "fort.31";
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
	
	if ( mrconf.Refs.length()==0 && !mrconf.autoRef )
	{
		cout << "Error: non optional parameter \"RefConfs\" missing." << endl;
		errors++;
	}
	
	if ( mrconf.Multiplicity == -1 )
	{
		cout << "Error: non optional parameter \"Multiplicity\" missing." << endl;
		errors++;
	}
	
	if ( mrconf.firstGuessConfs==-1 )
	{
		*outStream << "\"FirstGuessConfs\" is missing, defaults to 3." << endl;
		mrconf.firstGuessConfs = 3;
	}
	
	if ( mrconf.maxRefOpenShells==-1 )
	{
		*outStream << "\"maxRefOpenShells\" is missing, defaults to 4." << endl;
		mrconf.maxRefOpenShells = 4;
	}

	if ( mrconf.activeReferenceThreshold == -1 )
	{
		*outStream << "\"activeReferenceThreshold\" is missing, defaults to 0." << endl;
		mrconf.activeReferenceThreshold = 0;
	}
	
	if ( mrconf.nSelectionThresholds == 0 )
	{
		cout << "Error: non optional parameter \"SelectionThresholds\" missing." << endl;
		errors++;
	}

	if ( mrconf.NumberOfRootsE>0 && mrconf.NumberOfRoots != mrconf.NumberOfRootsE )
	{
		cout << "Error: number of roots in \"Roots\"- and \"RootEnergies\"-statement differ." << endl;
		errors++;
	}
	
	if ( mrconf.NumberOfRootsE>0 && 0 )
	{
		cout << "Error: in recalculation mode no specification of and \"RootEnergies\" allowed." << endl;
		errors++;
	}
	

FILE	*f;
	if (mrconf.MOIntegralFile.format != RIFormat)
	{
		if ( (f=fopen(mrconf.MOIntegralFile.name, "r"))==NULL )
		{
			cout << "Error: no file \"" << mrconf.MOIntegralFile.name
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



	qsort(mrconf.SelectionThreshold, 
		mrconf.nSelectionThresholds, sizeof(EnergyType),
		compThresh);

	mrconf.mrmos = new MRMOs(mrconf.MOIntegralFile);

	
INT	els = mrconf.NumberOfElectrons;
INT	elsErrors = 0;


	*outStream << endl;



	if ( mrconf.autoEquiv )
	{
//		if ( verbosity.isActive(Verbosity::DegenGuess) )
		{
			cout << "performing MO degeneration guess" << endl;
			cout << endl;
		}
		
	DiagIntContainer	diagInts(
		mrconf.getMOIntegralFile(), mrconf.getMultiplicity());

		mrconf.moEquivalence = diagInts.calcMODegeneration(
			mrconf.getNumberOfElectrons());

//		if ( verbosity.isActive(Verbosity::DegenGuess) )
		{
			cout << endl;
			cout << "degenerated MOs: " << *mrconf.moEquivalence << endl;
			cout << endl;
			cout << "MO degeneration guess finished" << endl;
			cout << endl;
			cout << endl;
		}
	}




	if ( mrconf.autoRef )
	{
		if ( verbosity.isActive(Verbosity::RefGuess) )
		{
			cout << "performing reference configuration guess" << endl;
			cout << endl;
		}
		
	DiagIntContainer	diagInts(
		mrconf.getMOIntegralFile(), mrconf.getMultiplicity());

	Configuration<MOType>	ground;
	
		mrconf.Refs = diagInts.calcRefConfs(
			mrconf.getNumberOfElectrons(),
			mrconf.getIrRep(),
			ground,
			mrconf.getNumberOfRoots()+mrconf.firstGuessConfs,
			mrconf.creatorSpace.length()>0);

//		mrconf.Refs.clear();
//		mrconf.Refs.add(ground);
		

		if ( mrconf.moRestriction )
			delete mrconf.moRestriction;
			
		mrconf.moRestriction = new MORestriction(
			ground, *mrconf.mrmos, 4);
			
		cout << *mrconf.moRestriction << endl;

		if ( verbosity.isActive(Verbosity::RefGuess) )
		{
			cout << "reference configuration guess finished" << endl;
			cout << endl;
			cout << endl;
		}
	}


// activate symmetry filter
//	if ( mrconf.creatorSpace.length()>0 )
	{
		mrconf.mrmos->initIntExt();
		mrconf.Refs = mrconf.genRefSpace(
			mrconf.Refs, mrconf.activeSpaceExcitationLevel);
	}

	if ( mrconf.PTRefs.length()==0 )
	{
		*outStream << "\"PTRefConfs\" is missing, defaults to \"RefConfs\"-input." << endl;
		mrconf.PTRefs = mrconf.Refs;
	}
	




	if ( mrconf.moRestriction )
		mrconf.moRestriction->setMaxMO(mrconf.mrmos->getMaxMO());

	if ( mrconf.moStatistics )
		mrconf.moStatistics = new MOStatistics(mrconf.mrmos->getMaxMO());


	if ( mrconf.moEquivalence && !mrconf.creatorSpace.length() )
	{
		cout << "symmetrizing reference space according to MO equivalences..." << endl;
	INT	n = mrconf.Refs.length();
		mrconf.moEquivalence->symmetrize(mrconf.Refs, mrconf.mrmos);
		cout << mrconf.Refs.length()-n << " configurations added." << endl;
		cout << endl;

		cout << "symmetrizing PT reference space according to MO equivalences..." << endl;
		n = mrconf.PTRefs.length();
		mrconf.moEquivalence->symmetrize(mrconf.PTRefs, mrconf.mrmos);
		cout << mrconf.PTRefs.length()-n << " configurations added." << endl;
		cout << endl;
		cout << endl;
	}


	{
	ConfigurationSet	set;
	Pix	i = mrconf.Refs.first();
	
		while ( i )
		{
		INT	n = mrconf.Refs(i).getNumberOfOpenShells();
			if ( n>mrconf.maxRefOpenShells || n<mrconf.Multiplicity-1 )
				set.add(mrconf.Refs(i));
			mrconf.Refs.next(i);
		}
		mrconf.Refs -= set;
	}
	
	{
	ConfigurationSet	set;
	Pix	i = mrconf.PTRefs.first();
	
		while ( i )
		{
		INT	n = mrconf.PTRefs(i).getNumberOfOpenShells();
			if ( n>mrconf.maxRefOpenShells || n<mrconf.Multiplicity-1 )
				set.add(mrconf.PTRefs(i));
			mrconf.PTRefs.next(i);
		}
		mrconf.PTRefs -= set;
	}
	
	mrconf.checkConfs(mrconf.Refs, els, elsErrors, errors);
	mrconf.checkConfs(mrconf.PTRefs, els, elsErrors, errors);


	mrconf.NumberOfElectrons = els;
	
	errors += elsErrors;


	mrconf.initInternalExternal();
	
	
	
	*outStream << endl;
	if ( errors )
	{
		cout << endl << errors << " errors.\naborting." << endl;
		exit(1);
	}
        delete outStream;
	return s;
}

void	MRConfInput::genActiveSpace(Configuration<MOType> conf, 
	INT *active, SLList<INT> down, SLList<INT> up)
{
INT	nIrreps = mrmos->getNumberOfIrReps();
INT	highest[nIrreps];


	for ( IrRep i=0 ; i<nIrreps ; i++ )
		highest[i] = mrmos->getStartMO(i);// - 1;
	
	for ( INT i=0 ; i<conf.getNumberOfOpenShells() ; i++ )
	{
		highest[mrmos->getIrRep(conf.getOpenShell(i))] =
			conf.getOpenShell(i);
		active[conf.getOpenShell(i)] = 1;
	}
	for ( INT i=0 ; i<conf.getNumberOfClosedShells() ; i++ )
		if ( highest[mrmos->getIrRep(conf.getClosedShell(i))] < 
					conf.getClosedShell(i) )
			highest[mrmos->getIrRep(conf.getClosedShell(i))] =
				conf.getClosedShell(i);
	
	
Pix	ai = down.first();		
Pix	ci = up.first();		
	for ( IrRep i=0 ; i<nIrreps ; i++ )
	{
		if ( ai )
		{
		INT	activeOccupied = down(ai);
			for ( MOType j=highest[i] ; 
				j>highest[i]-activeOccupied && j>0 && mrmos->getIrRep(j)==i ; j-- )
				active[j] = 1;
		}

		if ( ci )
		{
		INT	activeVirtual = up(ci);
			for ( MOType j = highest[i] ;
				j<highest[i]+activeVirtual && j<=mrmos->getMaxMO() && mrmos->getIrRep(j)==i ; j++ )
				active[j] = 1;
		}
			
		if ( ai )
			down.next(ai);		
		if ( ci )
			up.next(ci);		
	}
}



class Excitations {
public:
	Excitations(
		const MRMOs *mrmos, 
		INT level, 
		IrRep irrep,
		INT minOpenShells,
		Configuration<MOType> base,
		const INT *annihilatorMOs, const INT *creatorMOs);

	void	performExcitation(
		INT toAnnihilate, INT toCreate, 
		MOType mo1, MOType mo2, 
		Configuration<MOType> conf);

	ConfigurationSet	getConfSet() const {	return confs;	}

private:
const MRMOs *mrmos;
INT	level;
IrRep	irrep;
INT minOpenShells;
Configuration<MOType> base;
const INT *annihilatorMOs;
const INT *creatorMOs;
ConfigurationSet	confs;
INT	maxMO;
INT nIrreps;
};



Excitations::Excitations(
	const MRMOs *mrmos, 
	INT level,
	IrRep irrep,
	INT minOpenShells,
	Configuration<MOType> base,
	const INT *annihilatorMOs, const INT *creatorMOs) :
	mrmos(mrmos), level(level), irrep(irrep), minOpenShells(minOpenShells),
	base(base),
	annihilatorMOs(annihilatorMOs), creatorMOs(creatorMOs)
{
	maxMO = mrmos->getMaxMO();
	nIrreps = mrmos->getNumberOfIrReps();
	performExcitation(level, level, 1, 1, base);
}



void	Excitations::performExcitation(
	INT toAnnihilate, INT toCreate, 
	MOType mo1, MOType mo2, 
	Configuration<MOType> conf)
{
	if ( toAnnihilate )
	{
		for ( INT i=0 ; i<conf.getNumberOfOpenShells() ; i++ )
		{
		MOType	mo = conf.getOpenShell(i);

			if ( !annihilatorMOs[mo] || mo<mo1 )
				continue;

			conf.annihilate(mo);
			performExcitation(toAnnihilate-1, toCreate, mo, 1, conf);
			conf.create(mo);
		}
		for ( INT i=0 ; i<conf.getNumberOfClosedShells() ; i++ )
		{
		MOType	mo = conf.getClosedShell(i);

			if ( !annihilatorMOs[mo] || mo<mo1 )
				continue;

			conf.annihilate(mo);
			performExcitation(toAnnihilate-1, toCreate, mo, 1, conf);
			conf.create(mo);
		}
		return;
	}
	if ( toCreate )
	{
	MOType	mos = mo2;
	MOType	moe = maxMO;
	
		if ( toCreate==1 )
		{
		IrRep	matchIrrep = mrmos->getProd(
				conf.calcIrRep(*mrmos), irrep);

			if ( matchIrrep>=mrmos->getNumberOfIrReps() )
				return;
//			cout << "-----------" << endl;
//			cout << conf << endl;		
//			cout << "A: " << conf.calcIrRep(*mrmos) << "*" << irrep << " " << matchIrrep  << endl;

			mos = mrmos->getStartMO(matchIrrep);
			if ( matchIrrep+1<nIrreps )
				moe = mrmos->getStartMO(matchIrrep+1) - 1;
			else
				moe = maxMO;

			if ( mos<mo2 )
				mos = mo2;
//			cout << "mos=" << mos << ", moe=" << moe << endl;
		}
		for ( MOType mo=mos ; mo<=moe ; mo++ )
		{
			if ( creatorMOs[mo] )
			{
				if ( !conf.create(mo) )
				{
					performExcitation(0, toCreate-1, 0, mo, conf);
					conf.annihilate(mo);
				}
			}
		}
		return;
	}

	if ( conf.getNumberOfOpenShells()<minOpenShells )
		return;
//	cout << ":::::" << conf << endl;
	confs.add(conf);
}




ConfigurationSet	MRConfInput::genRefSpace(ConfigurationSet references, INT maxExc)
{
ConfigurationSet	newSet;
AVLSet<Configuration<MOType> >* internal = 		// internal configuration rests
	new AVLSet<Configuration<MOType> >[maxExc+1];						// with start addresses

INT	a[mrmos->getMaxMO()+1];
	memset(a, 0, (mrmos->getMaxMO()+1)*sizeof(INT));
INT	c[mrmos->getMaxMO()+1];
	memset(c, 0, (mrmos->getMaxMO()+1)*sizeof(INT));

/*
SLList<INT>	null;
Pix	iRef = references.first();
	while ( iRef )
	{
	Configuration<MOType>	conf(references(iRef));
		genActiveSpace(conf, a, annihilatorSpace, null);
		genActiveSpace(conf, c, null, creatorSpace);
		references.next(iRef);
	}
*/
/*	for ( INT i=1 ; i<=mrmos->getMaxMO() ; i++ )
	{
		if ( a[i] )
			mrmos->setInternal(i);
		if ( c[i] )
			mrmos->setInternal(i);
	}
	mrmos->initIntExt();
	
	cout << *mrmos << endl;
*/


	{
	Pix	i = annihilatorSpace.first();
		while ( i )
		{
			a[annihilatorSpace(i)] = 1;
			annihilatorSpace.next(i);
		}
	}
	{
	Pix	i = creatorSpace.first();
		while ( i )
		{
			c[creatorSpace(i)] = 1;
			creatorSpace.next(i);
		}
	}
	
	
	
	cout << endl << "active annihilators:" << endl;
	for ( INT i=1 ; i<=mrmos->getMaxMO() ; i++ )
		if ( a[i] )
			cout << i << " ";
	cout << endl;	

	cout << endl << "active creators:" << endl;
	for ( INT i=1 ; i<=mrmos->getMaxMO() ; i++ )
		if ( c[i] )
			cout << i << " ";
	cout << endl;	


Pix	iRef = references.first();
	while ( iRef )
	{
		if ( references(iRef).calcIrRep(*mrmos)==irrep )
			newSet.add(references(iRef));

		for ( INT i=1 ; i<=maxExc ; i++ )
		{
		Excitations	exc(mrmos, i, irrep, Multiplicity-1, references(iRef), a, c);
		ConfigurationSet	set(exc.getConfSet());
//		cout << "i=" << i << endl;
//		cout << set << endl << endl << endl;
		Pix	j = set.first();
			while ( j )
			{
				newSet.add(set(j));
				set.next(j);
			}
//			cout << newSet << endl;
		}
		references.next(iRef);
	}

/*	// use annihilators
Pix	iRef = references.first();
	while ( iRef )
	{
//			cout << "R: " << references(iRef) << endl;
		for ( INT i=0 ; i<=maxExc ; i++ )
		{
		MOIterator	Annihilators(i, mrmos, 1);
			while ( !Annihilators.isEnd() )
			{
			Configuration<MOType>	conf(references(iRef));
			INT	flag = 1;
			
				for ( INT j=0 ; j<i ; j++ )
					flag &= !!a[Annihilators.getMO(j)];
			
//				cout << ":::" << Annihilators << " " << flag << endl;
				if ( flag )
				{
					conf -= Annihilators;
					if ( conf.getNumberOfElectrons() ) //&& Annihilators.getNumberOfOpenShells()==0 )
						if ( i || conf.calcIrRep(*mrmos)==irrep )
							internal[i].add(conf);
				}
				Annihilators.next();
			}
		}
		references.next(iRef);
	}

	for ( INT i=0 ; i<=maxExc ; i++ )
	{
		cout << "::::::::::::::::: i=" << i << endl;
	Pix	iInt = internal[i].first();
		while ( iInt )
		{
			cout << ((Configuration<MOType>) internal[i](iInt)) << endl;
			internal[i].next(iInt);
		}
	}
	cout << "-------------------------------------" << endl;

	// use creators
	for ( INT i=0 ; i<=maxExc ; i++ )
	{
		cout << "********* i= " << i << endl;
	Pix	iInt = internal[i].first();
		while ( iInt )
		{
//		Configuration<MOType>	con(internal[i](iInt));
//			cout << "R: " << con << endl;
			MOIterator	Creators(i, mrmos, 1,
					mrmos->getProd(internal[i](iInt).calcIrRep(*mrmos), irrep));

			while ( !Creators.isEnd() )
			{
			Configuration<MOType>	conf(internal[i](iInt));
			
			INT	flag = 1;
			
				for ( INT j=0 ; j<i ; j++ )
					flag &= !!c[Creators.getMO(j)];
			
//				cout << ":::" << Creators << " " << flag << endl;
				if ( flag )
				{
					cout << conf << endl;
					cout << Creators << endl;
					conf += Creators;
					cout << conf << endl << endl << endl << "----------" << endl << endl;
					if ( conf.getNumberOfElectrons() ) //&& Creators.getNumberOfOpenShells()==0)
						newSet.add(conf);
				}
				Creators.next();
			}
			internal[i].next(iInt);
		}
	}

	for ( INT i=0 ; i<=maxExc ; i++ )
	{
		cout << "::::::::::::::::: i=" << i << endl;
	Pix	iInt = internal[i].first();
		while ( iInt )
		{
			cout << ((Configuration<MOType>) internal[i](iInt)) << endl;
			internal[i].next(iInt);
		}
	}
*/
        delete[] internal;
	return newSet;
}


		
/*
void	MRConfInput::setRefs.confs(const ConfigurationSet &confSet)
{
	if ( Refs.confs )
	{
		for ( INT i=0 ; i<Refs.n ; i++ )
			delete Refs.confs[i];
		free(Refs.confs);
	}
	
	
	Refs.n = confSet.length();
	
	Refs.confs = (Configuration<MOType> **) 
		malloc(Refs.n * sizeof(Configuration<MOType> *));

Pix	i = confSet.first(); 
INT	j = 0;
	while ( i )
	{
		Refs.confs[j++] = new Configuration<MOType>(*confSet(i));
		confSet.next(i);
	}

}
*/

void	MRConfInput::setVerbosity(const char *yytext)
{
INT	i = strlen(yytext) - 2;
	while ( i>=0 && yytext[i]!='=' )
		i--;
	i++;

istringstream	s(yytext+i, istringstream::in | istringstream::out);
	verbosity = Verbosity(s);
}


void	MRConfInput::setMORestriction(const char *yytext)
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
	moRestriction = new MORestriction(String(yytext+i));
}


void	MRConfInput::setMOEquivalence(const char *yytext, INT auto1)
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
