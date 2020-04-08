%{
#include <stdio.h>
#include <sstream>
#include "../../../lib/Container/String.h"
#include <ctype.h>
#include "RefSelInput.h"
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
istringstream	str(s,istringstream::in | istringstream::out);\
	str.seekg(0);\
	str >> conf;\
	ConfSet.add(conf);\
}

#define	ScanRoot(x) \
{\
istringstream	str(yytext,istringstream::in | istringstream::out);\
	str.seekg(0);\
	(x)->rootNumbers = (INT *) realloc((x)->rootNumbers, \
		++(x)->NumberOfRoots * sizeof(INT));\
	str >> (x)->rootNumbers[(x)->NumberOfRoots-1];\
}

#define	ScanRootEnergy(x) \
{\
istringstream	str(yytext,istringstream::in | istringstream::out);\
	str.seekg(0);\
	(x)->rootEnergies = (EnergyType *) realloc((x)->rootEnergies, \
		++(x)->NumberOfRootsE * sizeof(EnergyType));\
	str >> (x)->rootEnergies[(x)->NumberOfRootsE-1];\
}

#define	ScanThresh(x) \
{\
istringstream	str(yytext,istringstream::in | istringstream::out);\
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
istringstream	str(yytext,istringstream::in | istringstream::out);\
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
%s ScanFormat

	extern RefSelInput	*LexMRConf;
	extern RefSelInput	*LexPTMRConf;
	extern INT	LexMRConfErrors;


	SLList<INT>	*intList;
%%
.*#.*\n { INT i=strlen(yytext)-1; while ( i>=0 && yytext[i]!='#' ) i--; unput('\n'); for ( i-- ; i>=0 ; i-- ) unput(yytext[i]); }
<*>{ws}+	/* skip blanks and tabs */
<*>{separator}+		/* ignore blank separators */

MOIntegralFilename{ws}*={ws}*{name}{ws}*{separator}			ScanString(LexMRConf->MOIntegralFile.name);
MOIntegralFileFormat{ws}*=									BEGIN(ScanFormat);
FirstGuessConfs{ws}*={ws}*{posnum}{ws}*{separator}			ScanInt(LexMRConf->firstGuessConfs);
Verbosity{ws}*=.*											LexMRConf->setVerbosity(yytext);
ReferenceThreshold{ws}*={ws}*{floatnum}{ws}*{separator}				ScanDouble(LexMRConf->ReferenceThreshold);
MOEquivalence{ws}*={ws}*none{ws}*{separator}						LexMRConf->setMOEquivalence(NULL);
MOEquivalence{ws}*={ws}*auto{ws}*{separator}							LexMRConf->autoEquiv = 1;
MOEquivalence{ws}*=.*										LexMRConf->setMOEquivalence(yytext);
RefSelMode{ws}*={ws}*ConfThresh{ws}*{separator}						LexMRConf->refSelMode = RefSelMode::ConfThresh;
RefSelMode{ws}*={ws}*SumThresh{ws}*{separator}						LexMRConf->refSelMode = RefSelMode::SumThresh;
IrReps{ws}*={ws}*											BEGIN(IntSetPreScan); intList = &LexMRConf->irreps;



<RootScan,IntSetScan>\}	BEGIN(INITIAL);

<IntSetPreScan>\{{ws}*										BEGIN(IntSetScan);
<IntSetPreScan>.												cout << "unexpected token '" << yytext << "', expected '{'" << endl; LexMRConfErrors++;

<IntSetScan>[0-9]+												ScanIntSet;
<IntSetScan>[^ \t;\n]+											cout << "unexpected number '" << yytext << "' while scanning integer" << endl; LexMRConfErrors++;


<ScanFormat>{ws}*Old{ws}*{separator}						LexMRConf->MOIntegralFile.format=Fort31RecordFormatOld; BEGIN(INITIAL);
<ScanFormat>{ws}*New{ws}*{separator}						LexMRConf->MOIntegralFile.format=Fort31RecordFormatNew; BEGIN(INITIAL);
<ScanFormat>{ws}*TRADPT{ws}*{separator}						LexMRConf->MOIntegralFile.format=Fort31RecordFormatTRADPT; BEGIN(INITIAL);
<ScanFormat>{ws}*auto{ws}*{separator}						LexMRConf->MOIntegralFile.format=Fort31RecordFormatAuto; BEGIN(INITIAL);
<ScanFormat>{ws}*RIFormat{ws}*{separator}					LexMRConf->MOIntegralFile.format=RIFormat; BEGIN(INITIAL);
<ScanFormat>{ws}*GAMESSC1{ws}*{separator}					LexMRConf->MOIntegralFile.format=GAMESSC1; BEGIN(INITIAL);


<*>{TOKEN}													cout << "syntax error: unexpected token '" << yytext << "'" << endl; LexMRConfErrors++;

%%

