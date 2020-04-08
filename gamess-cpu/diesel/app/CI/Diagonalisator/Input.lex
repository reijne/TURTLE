%option noyywrap
%{
#include <stdio.h>
#include <sstream>
#include "../../../lib/Container/String.h"
#include <ctype.h>
#include "DiagInput.h"
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


#define	ScanRoot(x) \
{\
istringstream	str(yytext,istringstream::in | istringstream::out);\
	str.seekg(0);\
	(x)->rootNumbers = (INT *) realloc((x)->rootNumbers, \
		++(x)->NumberOfRoots * sizeof(INT));\
	str >> (x)->rootNumbers[(x)->NumberOfRoots-1];\
}


%}

ws	[ \t]

string		[ -}]*
alpha		[A-Za-z./]
dig			[0-9]
name		{alpha}({alpha}|{dig})+
posnum		{dig}+
floatnum	[-+]?{dig}*(\.{dig}+)?([eE][-+]?{dig}+)?
bool		"yes"|"no"

separator	[;\n]
TOKEN	[^, \n]+

%s RootPreScan
%s RootScan
%s ScanFormat

	extern DiagInput	*LexDiag;
	extern INT	LexDiagErrors;
%%
.*#.*\n { INT i=strlen(yytext)-1; while ( i>=0 && yytext[i]!='#' ) i--; unput('\n'); for ( i-- ; i>=0 ; i-- ) unput(yytext[i]); }
<*>{ws}+	/* skip blanks and tabs */
<*>{separator}+		/* ignore blank separators */

TempDir{ws}*={ws}*{string}{ws}*{separator}							ScanString(TempDir);
MOIntegralFilename{ws}*={ws}*{name}{ws}*{separator}					ScanString(LexDiag->MOIntegralFile.name);
MOIntegralFileFormat{ws}*=											BEGIN(ScanFormat);
ReferenceThreshold{ws}*={ws}*{floatnum}{ws}*{separator}				ScanDouble(LexDiag->RefThreshold);
MaxIters{ws}*={ws}*{posnum}{ws}*{separator}							ScanInt(LexDiag->MaxIters);
MaxHamiltonStorageMem{ws}*={ws}*{posnum}{ws}*GB{ws}*{separator}		ScanInt(LexDiag->MaxStorageMem); LexDiag->MaxStorageMem <<= 30;
MaxHamiltonStorageMem{ws}*={ws}*{posnum}{ws}*MB{ws}*{separator}		ScanInt(LexDiag->MaxStorageMem); LexDiag->MaxStorageMem <<= 20;
MaxHamiltonStorageMem{ws}*={ws}*{posnum}{ws}*KB{ws}*{separator}		ScanInt(LexDiag->MaxStorageMem); LexDiag->MaxStorageMem <<= 10;
MaxHamiltonStorageMem{ws}*={ws}*{posnum}{ws}*{separator}			ScanInt(LexDiag->MaxStorageMem);
PTReferenceThreshold{ws}*={ws}*{floatnum}{ws}*{separator}			ScanDouble(LexDiag->PTRefThreshold);
ConvergenceEnergyChange{ws}*={ws}*{floatnum}{ws}*{separator}		ScanDouble(LexDiag->ConvergenceEnergyChange);
ConvergenceEigenvectorChange{ws}*={ws}*{floatnum}{ws}*{separator}	ScanDouble(LexDiag->ConvergenceEigenvectorChange);
ConfTreeFileName{ws}*={ws}*{string}{ws}*{separator}                 ScanString(LexDiag->ConfTreeFileName);
StorePTEnergy{ws}*={ws}*no{ws}*{separator}							LexDiag->storePTEnergy = 0;
StorePTEnergy{ws}*={ws}*yes{ws}*{separator}							LexDiag->storePTEnergy = 1;
RootHoming{ws}*={ws}*no{ws}*{separator}								LexDiag->rootHoming = 0;
RootHoming{ws}*={ws}*yes{ws}*{separator}							LexDiag->rootHoming = 1;
StorePTCoef{ws}*={ws}*no{ws}*{separator}							LexDiag->storePTCoef = 0;
StorePTCoef{ws}*={ws}*yes{ws}*{separator}							LexDiag->storePTCoef = 1;
IterationMode{ws}*={ws}*CI{ws}*{separator}							LexDiag->iterationMode = DavidsonCI<double, double>::CI;
IterationMode{ws}*={ws}*ACPF{ws}*{separator}						LexDiag->iterationMode = DavidsonCI<double, double>::ACPF;
IterationMode{ws}*={ws}*AQCC{ws}*{separator}						LexDiag->iterationMode = DavidsonCI<double, double>::AQCC;
Precision{ws}*={ws}*float{ws}*{separator}							LexDiag->precision = DiagInput::floatPrec;
Precision{ws}*={ws}*double{ws}*{separator}							LexDiag->precision = DiagInput::doublePrec;
Verbosity{ws}*=.*													LexDiag->setVerbosity(yytext);

Roots{ws}*={ws}*													BEGIN(RootPreScan);

<RootScan>\}														BEGIN(INITIAL);

<RootPreScan>\{{ws}*												BEGIN(RootScan);
<RootPreScan>.														cout << "unexpected token '" << yytext << "', expected '{'" << endl; LexDiagErrors++;

<RootScan>[0-9]+													ScanRoot(LexDiag);
<RootScan>[^ \t;\n]+												cout << "unexpected number '" << yytext << "' while scanning roots" << endl; LexDiagErrors++;


<ScanFormat>{ws}*Old{ws}*{separator}								LexDiag->MOIntegralFile.format=Fort31RecordFormatOld; BEGIN(INITIAL);
<ScanFormat>{ws}*New{ws}*{separator}								LexDiag->MOIntegralFile.format=Fort31RecordFormatNew; BEGIN(INITIAL);
<ScanFormat>{ws}*TRADPT{ws}*{separator}								LexDiag->MOIntegralFile.format=Fort31RecordFormatTRADPT; BEGIN(INITIAL);
<ScanFormat>{ws}*auto{ws}*{separator}								LexDiag->MOIntegralFile.format=Fort31RecordFormatAuto; BEGIN(INITIAL);
<ScanFormat>{ws}*RIFormat{ws}*{separator}							LexDiag->MOIntegralFile.format=RIFormat; BEGIN(INITIAL);
<ScanFormat>{ws}*GAMESSC1{ws}*{separator}							LexDiag->MOIntegralFile.format=GAMESSC1; BEGIN(INITIAL);


<*>{TOKEN}															cout << "syntax error: unexpected token '" << yytext << "'" << endl; LexDiagErrors++;

%%
