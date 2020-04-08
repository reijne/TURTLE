%{
#include <stdio.h>
#include <sstream>
#include "../../../lib/Container/String.h"
#include <ctype.h>
#include "MRConfInput.h"
#include "../../../lib/QM/IO/Fortran/Fort31File.h"

using namespace std;

	INT	i;
	INT	nMRConf, nnMRConf;
	char	*readln();

#define isSeparator(x) ((x)=='\n' || (x)==';')



#define	ScanString(x) \
{\
INT	i = strlen(yytext) - 2;\
	while ( i>=0 && yytext[i]!='=' )\
		i--;\
	i++;\
	while ( yytext[i]==' ' )\
		i++;\
	yytext[strlen(yytext)-1] = 0;\
	(x) = String(yytext+i);\
}

#define	ScanInt(x) \
{\
int 	x_2;\
INT	i = strlen(yytext) - 2;\
	while ( i>=0 && yytext[i]!='=' )\
		i--;\
	i++;\
	sscanf(yytext+i, "%d", &(x_2));\
	x=x_2;\
}

#define	ScanDouble(x) \
{\
INT	i = strlen(yytext) - 2;\
	while ( i>=0 && yytext[i]!='=' )\
		i--;\
	i++;\
	sscanf(yytext+i, "%lf", &(x));\
}


#define	ScanConf(ConfSet) \
{\
Configuration<MOType>	conf;\
String	s(yytext);\
istringstream	str(yytext, istringstream::in | istringstream::out);\
	str.seekg(0);\
	str >> conf;\
	ConfSet.add(conf);\
}

#define	ScanRoot(x) \
{\
istringstream	str(yytext, istringstream::in | istringstream::out);\
	str.seekg(0);\
	(x)->rootNumbers = (INT *) realloc((x)->rootNumbers, \
		++(x)->NumberOfRoots * sizeof(INT));\
	str >> (x)->rootNumbers[(x)->NumberOfRoots-1];\
}

#define	ScanRootEnergy(x) \
{\
istringstream	str(yytext, istringstream::in | istringstream::out);\
	str.seekg(0);\
	(x)->rootEnergies = (EnergyType *) realloc((x)->rootEnergies, \
		++(x)->NumberOfRootsE * sizeof(EnergyType));\
	str >> (x)->rootEnergies[(x)->NumberOfRootsE-1];\
}

#define	ScanThresh(x) \
{\
istringstream	str(yytext, istringstream::in | istringstream::out);\
	str.seekg(0);\
	(x)->SelectionThreshold = (EnergyType *) realloc((x)->SelectionThreshold, \
		++(x)->nSelectionThresholds * sizeof(EnergyType));\
	str >> (x)->SelectionThreshold[(x)->nSelectionThresholds-1];\
	(x)->SelectionThresholdString = (char **) realloc((x)->SelectionThresholdString, \
		(x)->nSelectionThresholds * sizeof(char *));\
	(x)->SelectionThresholdString[(x)->nSelectionThresholds-1] = new char[strlen(yytext)+1];\
	strcpy((x)->SelectionThresholdString[(x)->nSelectionThresholds-1], yytext);\
}

#define	ScanNthExcite(x) \
{\
istringstream	str(yytext, istringstream::in | istringstream::out);\
	str.seekg(0);\
INT	i;\
	str >> i;\
	(x)->selectNthExcitation |= (1<<(i-1));\
}

#define	ScanIntSet \
{\
INT h = atoi(yytext);\
	intList->append(h);\
}

%}

%option noyywrap

ws	[ \t]

alpha		[A-Za-z./]
dig			[0-9]
name		{alpha}({alpha}|{dig})+
posnum		{dig}+
floatnum	[-+]?{dig}*(\.{dig}+)?([eE][-+]?{dig}+)?
bool		"yes"|"no"

separator	[;\n]
TOKEN	[^, \n]+

%s IntSetPreScan
%s IntSetScan
%s ConfPreScan
%s ConfScan
%s PTConfPreScan
%s PTConfScan
%s RootPreScan
%s RootEnergyPreScan
%s RootScan
%s RootEnergyScan
%s NthExcitePreScan
%s NthExciteScan
%s ThreshPreScan
%s ThreshScan
%s ScanFormat

	extern MRConfInput	*LexMRConf;
	extern MRConfInput	*LexPTMRConf;
	extern INT	LexMRConfErrors;


	SLList<INT>	*intList;
%%
.*#.*\n { INT i=strlen(yytext)-1; while ( i>=0 && yytext[i]!='#' ) i--; unput('\n'); for ( i-- ; i>=0 ; i-- ) unput(yytext[i]); }
<*>{ws}+	/* skip blanks and tabs */
<*>{separator}+		/* ignore blank separators */

MOIntegralFilename{ws}*={ws}*{name}{ws}*{separator}			ScanString(LexMRConf->MOIntegralFile.name);
MOIntegralFileFormat{ws}*=									BEGIN(ScanFormat);
NumberOfElectrons{ws}*={ws}*{posnum}{ws}*{separator}		ScanInt(LexMRConf->NumberOfElectrons);
ExcitationLevel{ws}*={ws}*{posnum}{ws}*{separator}			ScanInt(LexMRConf->ExcitationLevel);
Multiplicity{ws}*={ws}*{posnum}{ws}*{separator}				ScanInt(LexMRConf->Multiplicity);
IrRep{ws}*={ws}*{posnum}{ws}*{separator}					ScanInt(LexMRConf->irrep);
ActiveSpaceExcitationLevel{ws}*={ws}*{posnum}{ws}*{separator}					ScanInt(LexMRConf->activeSpaceExcitationLevel);
maxRefOpenShells{ws}*={ws}*{posnum}{ws}*{separator}					ScanInt(LexMRConf->maxRefOpenShells);
SelectInternal{ws}*={ws}*no{ws}*{separator}					LexMRConf->selectInternal = 0;
SelectInternal{ws}*={ws}*yes{ws}*{separator}				LexMRConf->selectInternal = 1;
StorePTEnergy{ws}*={ws}*no{ws}*{separator}					LexMRConf->storePTEnergy = 0;
StorePTEnergy{ws}*={ws}*yes{ws}*{separator}					LexMRConf->storePTEnergy = 1;
StorePTCoef{ws}*={ws}*no{ws}*{separator}					LexMRConf->storePTCoef = 0;
StorePTCoef{ws}*={ws}*yes{ws}*{separator}					LexMRConf->storePTCoef = 1;
MOStatistics{ws}*={ws}*no{ws}*{separator}					LexMRConf->moStatistics = (MOStatistics *) 0;
MOStatistics{ws}*={ws}*yes{ws}*{separator}					LexMRConf->moStatistics = (MOStatistics *) 1;
RefConfs{ws}*={ws}*auto{ws}*{separator}						LexMRConf->autoRef = 1;
AnnihilatorSpace{ws}*={ws}*									BEGIN(IntSetPreScan); intList = &LexMRConf->annihilatorSpace;
CreatorSpace{ws}*={ws}*										BEGIN(IntSetPreScan); intList = &LexMRConf->creatorSpace;
RefConfs{ws}*={ws}*											BEGIN(ConfPreScan);
PTRefConfs{ws}*={ws}*										BEGIN(PTConfPreScan);
Roots{ws}*={ws}*											BEGIN(RootPreScan);
RootEnergies{ws}*={ws}*										BEGIN(RootEnergyPreScan);
SelectionThresholds{ws}*={ws}*								BEGIN(ThreshPreScan);
selectNthExcitation{ws}*={ws}*								BEGIN(NthExcitePreScan);
Verbosity{ws}*=.*											LexMRConf->setVerbosity(yytext);
MORestrictions{ws}*={ws}*none{ws}*							LexMRConf->setMORestriction(NULL);
MORestrictions{ws}*=.*										LexMRConf->setMORestriction(yytext);
MOEquivalence{ws}*={ws}*none{ws}*							LexMRConf->setMOEquivalence(NULL);
MOEquivalence{ws}*={ws}*auto{ws}*							LexMRConf->autoEquiv = 1;
MOEquivalence{ws}*=.*										LexMRConf->setMOEquivalence(yytext);
EstimationMode{ws}*={ws}*EpsteinNesbet{ws}*{separator}		LexMRConf->estimationMode = EnergyMap::EpsteinNesbet;
EstimationMode{ws}*={ws}*Wenzel{ws}*{separator}				LexMRConf->estimationMode = EnergyMap::Wenzel;
ActiveReferenceThreshold{ws}*={ws}*{floatnum}{ws}*{separator}				ScanDouble(LexMRConf->activeReferenceThreshold);
RandomProb{ws}*={ws}*{floatnum}{ws}*{separator}				ScanDouble(LexMRConf->randomProb);


<ConfPreScan>\{{ws}*										BEGIN(ConfScan);
<ConfPreScan>.												cout << "unexpected token '" << yytext << "', expected '{'" << endl; LexMRConfErrors++;

<ConfScan>[0-9\- \t]+{separator}							ScanConf(LexMRConf->getRefConfSet());
<ConfScan>[^ 0-9\-\n\}]*									cout << "unexpected number '" << yytext << "' while scanning configurations" << endl; LexMRConfErrors++;


<PTConfPreScan>\{{ws}*										BEGIN(PTConfScan);
<PTConfPreScan>.											cout << "unexpected token '" << yytext << "', expected '{'" << endl; LexMRConfErrors++;

<PTConfScan>[0-9\- \t]+{separator}							ScanConf(LexMRConf->getPTRefConfSet());
<PTConfScan>[^ 0-9\-\n\}]*									cout << "unexpected number '" << yytext << "' while scanning configurations" << endl; LexMRConfErrors++;

<ConfScan,PTConfScan,RootScan,NthExciteScan,ThreshScan,RootEnergyScan,IntSetScan>\}	BEGIN(INITIAL);

<IntSetPreScan>\{{ws}*										BEGIN(IntSetScan);
<IntSetPreScan>.												cout << "unexpected token '" << yytext << "', expected '{'" << endl; LexMRConfErrors++;

<IntSetScan>[0-9]+												ScanIntSet;
<IntSetScan>[^ \t;\n]+											cout << "unexpected number '" << yytext << "' while scanning integer" << endl; LexMRConfErrors++;

<RootPreScan>\{{ws}*										BEGIN(RootScan);
<RootPreScan>.												cout << "unexpected token '" << yytext << "', expected '{'" << endl; LexMRConfErrors++;

<RootScan>[0-9]+											ScanRoot(LexMRConf);
<RootScan>[^ \t;\n]+										cout << "unexpected number '" << yytext << "' while scanning roots" << endl; LexMRConfErrors++;

<RootEnergyPreScan>\{{ws}*									BEGIN(RootEnergyScan);
<RootEnergyPreScan>.										cout << "unexpected token '" << yytext << "', expected '{'" << endl; LexMRConfErrors++;

<RootEnergyScan>{floatnum}									ScanRootEnergy(LexMRConf);
<RootEnergyScan>[^ \t;\n]+									cout << "unexpected number '" << yytext << "' while scanning root energies" << endl; LexMRConfErrors++;

<ThreshPreScan>\{{ws}*										BEGIN(ThreshScan);
<ThreshPreScan>.											cout << "unexpected token '" << yytext << "', expected '{'" << endl; LexMRConfErrors++;

<ThreshScan>{floatnum}										ScanThresh(LexMRConf);
<ThreshScan>[^ \t;\n]+										cout << "unexpected number '" << yytext << "' while scanning thresholds" << endl; LexMRConfErrors++;

<NthExcitePreScan>\{{ws}*									BEGIN(NthExciteScan);
<NthExcitePreScan>.											cout << "unexpected token '" << yytext << "', expected '{'" << endl; LexMRConfErrors++;

<NthExciteScan>[0-9]+										ScanNthExcite(LexMRConf);
<NthExciteScan>[^ \t;\n]+									cout << "unexpected number '" << yytext << "' while scanning selectNthExcitations" << endl; LexMRConfErrors++;




<ScanFormat>{ws}*Old{ws}*{separator}						LexMRConf->MOIntegralFile.format=Fort31RecordFormatOld; BEGIN(INITIAL);
<ScanFormat>{ws}*New{ws}*{separator}						LexMRConf->MOIntegralFile.format=Fort31RecordFormatNew; BEGIN(INITIAL);
<ScanFormat>{ws}*TRADPT{ws}*{separator}						LexMRConf->MOIntegralFile.format=Fort31RecordFormatTRADPT; BEGIN(INITIAL);
<ScanFormat>{ws}*auto{ws}*{separator}						LexMRConf->MOIntegralFile.format=Fort31RecordFormatAuto; BEGIN(INITIAL);
<ScanFormat>{ws}*RIFormat{ws}*{separator}					LexMRConf->MOIntegralFile.format=RIFormat; BEGIN(INITIAL);
<ScanFormat>{ws}*GAMESSC1{ws}*{separator}					LexMRConf->MOIntegralFile.format=GAMESSC1; BEGIN(INITIAL);


<*>{TOKEN}													cout << "syntax error: unexpected token '" << yytext << "'" << endl; LexMRConfErrors++;

%%

