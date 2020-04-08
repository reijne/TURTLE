%option noyywrap
%{
#include <stdio.h>
#include <sstream>
#include "../../../../lib/Container/String.h"
#include <ctype.h>
#include "DieselInput.h"
#include "../../../../lib/QM/IO/Fortran/Fort31File.h"

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
int	x_2;\
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


#define	ScanIntSet \
{\
INT h = atoi(yytext);\
	intList->append(h);\
}

#define	ScanIntVector \
{\
INT h = atoi(yytext);\
	intVector->push_back(h);\
}

#define	ScanDoubleSet \
{\
double h = atof(yytext);\
	doubleList->append(h);\
}

#define	ScanStringSet \
{\
String h(yytext);\
	stringList->append(h);\
}

#define	ScanConf(ConfSet) \
{\
Configuration<MOType>	conf;\
/*String	s(yytext); */\
istringstream	str(yytext,istringstream::in | istringstream::out);\
	str.seekg(0);\
	str >> conf;\
	ConfSet.add(conf);\
}

%}

ws	[ \t]

string		[ -}]*
alpha		[A-Za-z\.]
dig			[0-9]
name		{alpha}({alpha}|{dig})+
posnum		{dig}+
floatnum	[-+]?{dig}*(\.{dig}+)?([eE][-+]?{dig}+)?
bool		"yes"|"no"

separator	[;\n]
TOKEN	[^, \n]+

%s ConfPreScan
%s ConfScan
%s IntSetPreScan
%s IntSetScan
%s IntVectorPreScan
%s IntVectorScan
%s DoubleSetPreScan
%s DoubleSetScan
%s StringSetPreScan
%s StringSetScan

	extern DieselInput	*LexDiesel;
	extern INT	LexDieselErrors;
	
	SLList<INT>	*intList;
	vector<INT>	*intVector;
	SLList<double>	*doubleList;
	SLList<String>	*stringList;
%%
.*#.*\n { INT i=strlen(yytext)-1; while ( i>=0 && yytext[i]!='#' ) i--; unput('\n'); for ( i-- ; i>=0 ; i-- ) unput(yytext[i]); }
<*>{ws}+	/* skip blanks and tabs */
<*>{separator}+		/* ignore blank separators */

MOLCASRootDir{ws}*={ws}*{string}{ws}*{separator}			ScanString(LexDiesel->MOLCASRootDir);
TempDir{ws}*={ws}*{string}{ws}*{separator}					ScanString(LexDiesel->TempDir);
MOIntegralFilename{ws}*={ws}*{name}{ws}*{separator}			ScanString(LexDiesel->MOIntegralFile.name);
MOIntegralFileFormat{ws}*={ws}*Old{ws}*{separator}			LexDiesel->MOIntegralFile.format=Fort31RecordFormatOld;
MOIntegralFileFormat{ws}*={ws}*New{ws}*{separator}			LexDiesel->MOIntegralFile.format=Fort31RecordFormatNew;
MOIntegralFileFormat{ws}*={ws}*TRADPT{ws}*{separator}			LexDiesel->MOIntegralFile.format=Fort31RecordFormatTRADPT;
MOIntegralFileFormat{ws}*={ws}*auto{ws}*{separator}			LexDiesel->MOIntegralFile.format=Fort31RecordFormatAuto;
MOIntegralFileFormat{ws}*={ws}*RIFormat{ws}*{separator}			LexDiesel->MOIntegralFile.format=RIFormat;
MOIntegralFileFormat{ws}*={ws}*GAMESSC1{ws}*{separator}			LexDiesel->MOIntegralFile.format=GAMESSC1;


Roots{ws}*={ws}*											BEGIN(IntVectorPreScan); intVector = &LexDiesel->rootsDefault;
Roots0{ws}*={ws}*											BEGIN(IntVectorPreScan); intVector = &LexDiesel->roots[0];
Roots1{ws}*={ws}*											BEGIN(IntVectorPreScan); intVector = &LexDiesel->roots[1];
Roots2{ws}*={ws}*											BEGIN(IntVectorPreScan); intVector = &LexDiesel->roots[2];
Roots3{ws}*={ws}*											BEGIN(IntVectorPreScan); intVector = &LexDiesel->roots[3];
Roots4{ws}*={ws}*											BEGIN(IntVectorPreScan); intVector = &LexDiesel->roots[4];
Roots5{ws}*={ws}*											BEGIN(IntVectorPreScan); intVector = &LexDiesel->roots[5];
Roots6{ws}*={ws}*											BEGIN(IntVectorPreScan); intVector = &LexDiesel->roots[6];
Roots7{ws}*={ws}*											BEGIN(IntVectorPreScan); intVector = &LexDiesel->roots[7];
PreScanRoots{ws}*={ws}*{posnum}{ws}*{separator}	 		ScanInt(LexDiesel->preScanRoots);
SelectionThresholds{ws}*={ws}*								BEGIN(StringSetPreScan); stringList = &LexDiesel->thresholds;
MRPTSelectionThresholds{ws}*={ws}*							BEGIN(StringSetPreScan); stringList = &LexDiesel->MRPTthresholds;
NumberOfElectrons{ws}*={ws}*{posnum}{ws}*{separator}		ScanInt(LexDiesel->NumberOfElectrons);
ActiveReferenceThreshold{ws}*={ws}*{floatnum}{ws}*{separator}		ScanDouble(LexDiesel->activeReferenceThreshold);
ExcitationLevel{ws}*={ws}*{posnum}{ws}*{separator}			ScanInt(LexDiesel->ExcitationLevel);
maxRefGenIters{ws}*={ws}*{posnum}{ws}*{separator}			ScanInt(LexDiesel->maxRefGenIters);
ActiveSpaceExcitationLevel{ws}*={ws}*{posnum}{ws}*{separator}					ScanInt(LexDiesel->activeSpaceExcitationLevel);
maxRefOpenShells{ws}*={ws}*{posnum}{ws}*{separator}					ScanInt(LexDiesel->maxRefOpenShells);
Multiplicities{ws}*={ws}*									BEGIN(IntSetPreScan); intList = &LexDiesel->multiplicities;
IrReps{ws}*={ws}*											BEGIN(IntSetPreScan); intList = &LexDiesel->irreps;
SelectInternal{ws}*={ws}*no{ws}*{separator}					LexDiesel->selectInternal = 0;
SelectInternal{ws}*={ws}*yes{ws}*{separator}				LexDiesel->selectInternal = 1;
RefConfs{ws}*={ws}*auto{ws}*{separator}						LexDiesel->autoRef = 1;
RefConfs{ws}*={ws}*											BEGIN(ConfPreScan);
selectNthExcitation{ws}*={ws}*								BEGIN(IntSetPreScan); intList = &LexDiesel->nthExcitation;
AnnihilatorSpace{ws}*={ws}*									BEGIN(IntSetPreScan); intList = &LexDiesel->annihilatorSpace;
CreatorSpace{ws}*={ws}*										BEGIN(IntSetPreScan); intList = &LexDiesel->creatorSpace;
MORestrictions{ws}*={ws}*none{ws}*{separator}						LexDiesel->setMORestriction(NULL);
MORestrictions{ws}*=.*										LexDiesel->setMORestriction(yytext);
MOEquivalence{ws}*={ws}*none{ws}*{separator}						LexDiesel->setMOEquivalence(NULL);
MOEquivalence{ws}*={ws}*auto{ws}*{separator}							LexDiesel->autoEquiv = 1;
MOEquivalence{ws}*=.*										LexDiesel->setMOEquivalence(yytext);
SelectionEstimationMode{ws}*={ws}*EpsteinNesbet{ws}*{separator}		LexDiesel->estimationMode = EnergyMap::EpsteinNesbet;
SelectionEstimationMode{ws}*={ws}*Wenzel{ws}*{separator}				LexDiesel->estimationMode = EnergyMap::Wenzel;
MaxHamiltonStorageMem{ws}*={ws}*{posnum}{ws}*GB{ws}*{separator}		ScanInt(LexDiesel->MaxHamiltonStorageMem); LexDiesel->MaxHamiltonStorageMem <<= 30;
MaxHamiltonStorageMem{ws}*={ws}*{posnum}{ws}*MB{ws}*{separator}		ScanInt(LexDiesel->MaxHamiltonStorageMem); LexDiesel->MaxHamiltonStorageMem <<= 20;
MaxHamiltonStorageMem{ws}*={ws}*{posnum}{ws}*KB{ws}*{separator}		ScanInt(LexDiesel->MaxHamiltonStorageMem); LexDiesel->MaxHamiltonStorageMem <<= 10;
MaxHamiltonStorageMem{ws}*={ws}*{posnum}{ws}*{separator}			ScanInt(LexDiesel->MaxHamiltonStorageMem);
ConvergenceEnergyChange{ws}*={ws}*{floatnum}{ws}*{separator}		ScanDouble(LexDiesel->ConvergenceEnergyChange);
ConvergenceEigenvectorChange{ws}*={ws}*{floatnum}{ws}*{separator}	ScanDouble(LexDiesel->ConvergenceEigenvectorChange);
ReferenceThreshold{ws}*={ws}*{floatnum}{ws}*{separator}	ScanDouble(LexDiesel->RefThreshold);
PTReferenceThreshold{ws}*={ws}*{floatnum}{ws}*{separator}	ScanDouble(LexDiesel->PTRefThreshold);
NatOrbReferenceThreshold{ws}*={ws}*{floatnum}{ws}*{separator}	ScanDouble(LexDiesel->NatOrbRefThreshold);
MRPTInhomogenityThreshold{ws}*={ws}*{floatnum}{ws}*{separator}	ScanString(LexDiesel->MRPTInhomogenityThreshold);
nProcessors{ws}*={ws}*{posnum}{ws}*{separator}							ScanInt(LexDiesel->nProc);
MaxDavidsonIters{ws}*={ws}*{posnum}{ws}*{separator}							ScanInt(LexDiesel->MaxIters);
FirstGuessConfs{ws}*={ws}*{posnum}{ws}*{separator}							ScanInt(LexDiesel->firstGuessConfs);
useNaturalOrbitals{ws}*={ws}*no{ws}*{separator}					LexDiesel->useNaturalOrbitals = 0;
useNaturalOrbitals{ws}*={ws}*yes{ws}*{separator}					LexDiesel->useNaturalOrbitals = 1;
averagedNaturalOrbitals{ws}*={ws}*no{ws}*{separator}					LexDiesel->averagedNaturalOrbitals = 0;
averagedNaturalOrbitals{ws}*={ws}*yes{ws}*{separator}					LexDiesel->averagedNaturalOrbitals = 1;
NaturalOrbitalSelectionThreshold{ws}*={ws}*{floatnum}{ws}*{separator}		ScanString(LexDiesel->NaturalOrbitalSelectionThreshold);
fullMRCIExtrapolation{ws}*={ws}*											BEGIN(StringSetPreScan); stringList = &LexDiesel->fullMRCIExtrapolation;
propertyThresholds{ws}*={ws}*											BEGIN(StringSetPreScan); stringList = &LexDiesel->propertyThresholds;
orbitalFile{ws}*={ws}*{name}{ws}*{separator}			ScanString(LexDiesel->orbitalFile);
store{ws}*=.*											BEGIN(StringSetPreScan); stringList = &LexDiesel->store;
RefSelMode{ws}*={ws}*ConfThresh{ws}*{separator}						LexDiesel->refSelMode = RefSelMode::ConfThresh;
RefSelMode{ws}*={ws}*SumThresh{ws}*{separator}						LexDiesel->refSelMode = RefSelMode::SumThresh;

HomogeneousSelection{ws}*={ws}*no{ws}*{separator}					LexDiesel->homogeneousSelection = 0;
HomogeneousSelection{ws}*={ws}*yes{ws}*{separator}				LexDiesel->homogeneousSelection = 1;

Verbosity{ws}*=.*											LexDiesel->setVerbosity(yytext);
RootHoming{ws}*={ws}*no{ws}*{separator}								LexDiesel->rootHoming = 0;
RootHoming{ws}*={ws}*yes{ws}*{separator}							LexDiesel->rootHoming = 1;
UsePreviousSelection{ws}*={ws}*no{ws}*{separator}					LexDiesel->usePreviousSelection = 0;
UsePreviousSelection{ws}*={ws}*yes{ws}*{separator}					LexDiesel->usePreviousSelection = 1;



<ConfPreScan>\{{ws}*										BEGIN(ConfScan);
<ConfPreScan>.												cout << "unexpected token '" << yytext << "', expected '{'" << endl; LexDieselErrors++;

<ConfScan>[0-9\- \t]+{separator}							ScanConf(LexDiesel->getRefConfSet());
<ConfScan>[^ 0-9\-\n\}]*									cout << "unexpected number '" << yytext << "' while scanning configurations" << endl; LexDieselErrors++;


<ConfScan,IntSetScan,IntVectorScan,DoubleSetScan,StringSetScan>\}										BEGIN(INITIAL);

<IntSetPreScan>\{{ws}*										BEGIN(IntSetScan);
<IntSetPreScan>.												cout << "unexpected token '" << yytext << "', expected '{'" << endl; LexDieselErrors++;

<IntSetScan>[0-9]+												ScanIntSet;
<IntSetScan>[^ \t;\n]+											cout << "unexpected number '" << yytext << "' while scanning integer" << endl; LexDieselErrors++;

<IntVectorPreScan>\{{ws}*										BEGIN(IntVectorScan);
<IntVectorPreScan>.												cout << "unexpected token '" << yytext << "', expected '{'" << endl; LexDieselErrors++;

<IntVectorScan>[0-9]+												ScanIntVector;
<IntVectorScan>[^ \t;\n]+											cout << "unexpected number '" << yytext << "' while scanning integer" << endl; LexDieselErrors++;

<DoubleSetPreScan>\{{ws}*										BEGIN(DoubleSetScan);
<DoubleSetPreScan>.												cout << "unexpected token '" << yytext << "', expected '{'" << endl; LexDieselErrors++;

<DoubleSetScan>{floatnum}												ScanDoubleSet;
<DoubleSetScan>[^ \t;\n]+											cout << "unexpected number '" << yytext << "' while scanning double" << endl; LexDieselErrors++;

<StringSetPreScan>\{{ws}*										BEGIN(StringSetScan);
<StringSetPreScan>.												cout << "unexpected token '" << yytext << "', expected '{'" << endl; LexDieselErrors++;

<StringSetScan>[^ \t\}]*												ScanStringSet;
<StringSetScan>[^ \t;\n]+											cout << "unexpected token '" << yytext << "' while scanning string" << endl; LexDieselErrors++;



<*>{TOKEN}													cout << "syntax error: unexpected token '" << yytext << "'" << endl; LexDieselErrors++;

%%

