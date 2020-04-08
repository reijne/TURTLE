%{
#include <stdio.h>
#include <sstream>
#include "../../../lib/Container/String.h"
#include <ctype.h>
#include "MRPTInput.h"
#include "../../../lib/Container/TempDir.h"
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
INT	i = strlen(yytext) - 2;\
	while ( i>=0 && yytext[i]!='=' )\
		i--;\
	i++;\
	sscanf(yytext+i,FORMATI , &(x));\
}

#define	ScanDouble(x) \
{\
INT	i = strlen(yytext) - 2;\
	while ( i>=0 && yytext[i]!='=' )\
		i--;\
	i++;\
	sscanf(yytext+i, "%lf", &(x));\
}



#define	ScanRoot(x) \
{\
istringstream	str(yytext,istringstream::in | istringstream::out);\
	str.seekg(0);\
	(x)->rootNumbers = (INT *) realloc((x)->rootNumbers, \
		++(x)->NumberOfRoots * sizeof(INT));\
	str >> (x)->rootNumbers[(x)->NumberOfRoots-1];\
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


%}

%option noyywrap

ws	[ \t]

string		[ -}]*
alpha		[A-Za-z]
dig			[0-9]
name		{alpha}({alpha}|{dig})+
posnum		{dig}+
floatnum	[-+]?{dig}*(\.{dig}+)?([eE][-+]?{dig}+)?
bool		"yes"|"no"

separator	[;\n]
TOKEN	[^, \n]+

%s RootPreScan
%s RootScan
%s ThreshPreScan
%s ThreshScan
%s ScanFormat

	extern MRPTInput	*LexMRPT;
	extern INT	LexMRPTErrors;
%%
<*>{ws}+	/* skip blanks and tabs */
<*>^\#.*	/* ignore comments */
<*>{separator}+		/* ignore blank separators */

TempDir{ws}*={ws}*{string}{ws}*{separator}					ScanString(TempDir);
MOIntegralFilename{ws}*={ws}*{name}{ws}*{separator}			ScanString(LexMRPT->MOIntegralFile.name);
MOIntegralFileFormat{ws}*=									BEGIN(ScanFormat);
Roots{ws}*={ws}*											BEGIN(RootPreScan);
SelectionThresholds{ws}*={ws}*								BEGIN(ThreshPreScan);
ProjectionMode{ws}*={ws}*no0{ws}*{separator}				LexMRPT->projectionMode = MRMPH0Matrix<double, double>::P_no0;
ProjectionMode{ws}*={ws}*0Complement{ws}*{separator}		LexMRPT->projectionMode = MRMPH0Matrix<double, double>::P_0Complement;
calcMP3{ws}*={ws}*yes{ws}*{separator}						LexMRPT->calcMP3 = 1;
calcMP3{ws}*={ws}*no{ws}*{separator}						LexMRPT->calcMP3 = 0;
InhomogenityThreshold{ws}*={ws}*{name}{ws}*{separator}		ScanString(LexMRPT->inhomogenityThreshold);
InhomogenityThreshold{ws}*={ws}*{floatnum}{ws}*{separator}		ScanString(LexMRPT->inhomogenityThreshold);
Verbosity{ws}*=.*											LexMRPT->setVerbosity(yytext);

<RootScan,ThreshScan>\}										BEGIN(INITIAL);

<RootPreScan>\{{ws}*										BEGIN(RootScan);
<RootPreScan>.												cout << "unexpected token '" << yytext << "', expected '{'" << endl; LexMRPTErrors++;

<RootScan>[0-9]+											ScanRoot(LexMRPT);
<RootScan>[^ \t;\n]+										cout << "unexpected number '" << yytext << "' while scanning roots" << endl; LexMRPTErrors++;

<ThreshPreScan>\{{ws}*										BEGIN(ThreshScan);
<ThreshPreScan>.											cout << "unexpected token '" << yytext << "', expected '{'" << endl; LexMRPTErrors++;

<ThreshScan>{floatnum}										ScanThresh(LexMRPT);
<ThreshScan>[^ \t;\n]+										cout << "unexpected number '" << yytext << "' while scanning thresholds" << endl; LexMRPTErrors++;


<ScanFormat>{ws}*Old{ws}*{separator}						LexMRPT->MOIntegralFile.format=Fort31RecordFormatOld; BEGIN(INITIAL);
<ScanFormat>{ws}*New{ws}*{separator}						LexMRPT->MOIntegralFile.format=Fort31RecordFormatNew; BEGIN(INITIAL);
<ScanFormat>{ws}*TRADPT{ws}*{separator}						LexMRPT->MOIntegralFile.format=Fort31RecordFormatTRADPT; BEGIN(INITIAL);
<ScanFormat>{ws}*auto{ws}*{separator}						LexMRPT->MOIntegralFile.format=Fort31RecordFormatAuto; BEGIN(INITIAL);
<ScanFormat>{ws}*RIFormat{ws}*{separator}					LexMRPT->MOIntegralFile.format=RIFormat; BEGIN(INITIAL);
<ScanFormat>{ws}*GAMESSC1{ws}*{separator}					LexMRPT->MOIntegralFile.format=GAMESSC1; BEGIN(INITIAL);


<*>{TOKEN}													cout << "syntax error: unexpected token '" << yytext << "'" << endl; LexMRPTErrors++;

%%

