//***********************************************************************
//
//	Name:	TableCases.cc
//
//	Description:	implements configuration handling
//					based on second quantization
//
//	Author:	Michael Hanrath
//
//	Version:	0.0
//
//	Date:	08.06.1996
//
//	Literature:	Volker Pleﬂ: "Ein direktes, individuell
//			selektierendes Multireferenz Konfigurations-
//			Wechselwirkungsverfahren", 
//			Dissertation Uni Bonn, 1994
//
//
//
//***********************************************************************


#include "TableCase.h"


#include <string>
#include <stdlib.h>


static const INT	NumberOfCases = 29;	
//	Table 3.2 in Volker Pleﬂ Diss.
static const struct TCases {
	INT	dK;
	INT	P;
	INT	R;
	char *name;
	} Cases[NumberOfCases] = 
	{	
	// higher excitation
		{	0,	0,  0,	  "0"   },		//	0


	// double excitations
		{	0,	4,  1,    "1"	},		//	1

		{	0,	2,  1,    "2"	},		//	2 (!)
		{	0,	2,  1,    "2"	},		//	3 (!)
		{	1,	2,  2,    "3"	},		//	4
		{	1,	2,  1,    "4"	},		//	5

		{	0,	1,  1,    "5"	},		//	6
		{	0,	1,  6,    "6"	},		//	7

		{	0,	1,  2,   "7a"	},		//	8
		{	0,	1,  3,   "7b"	},		//	9
		{	0,	1,  4,   "7c"	},		//	10
		{	0,	1,  5,   "7d"	},		//	11

		{	1,	1,  4,   "8a"	},		//	12
		{	1,	1,  5,   "8b"	},		//	13
		{	1,	1,  6,   "8c"	},		//	14

		{	1,	1,  1,   "9a"	},		//	15
		{	1,	1,  2,   "9b"	},		//	16
		{	1,	1,  3,   "9c"	},		//	17

		{	2,	1,  1,  "10a"	},		//	18
		{	2,	1,  2,  "10b"	},		//	19
		{	2,	1,  3,  "10c"	},		//	20
		{	2,	1,  4,  "10d"	},		//	21
		{	2,	1,  5,  "10e"	},		//	22
		{	2,	1,  6,  "10f"	},		//	23


	// single excitations
		{	0,	3,  1,   "11"	},		//	24
		{	0,	3,  2,   "12"	},		//	25
		{	1,	3,  1,  "13a"	},		//	26
		{	1,	3,  2,  "13b"	},		//	27


	// no excitation
		{	0,	5,  1,   "14"	}		//	28
	};		


static unsigned char CoulombExchangeIntegral[NumberOfCases][4] =
	//
	//	binary coding of positions for coulomb- and exchange-integral:
	//
	//	 7		 6		 5		 4		 3		 2		 1		 0
	//	2		2		2		2		2		2		2		2
	// -------------------------------------------------------------
	//	1./2.	S/D		|--Pos.--| 		1./2.	S/D		|--Pos.--|
	//	|--------exchange--------|		|--------coulomb---------|
	//
	{	
        // higher excitation
                { 0x00, 0x00, 0x00, 0x00  },              //      0


        // double excitations
                { 0x04, 0x0c, 0x04, 0x0c  },              //      1

                { 0x04, 0x00, 0x04, 0x08  },              //      2 (!)
                { 0x0c, 0x04, 0x0c, 0x00  },              //      3 (!)
                { 0x04, 0x09, 0x04, 0x08  },              //      4
                { 0x0c, 0x09, 0x0c, 0x08  },              //      5

                { 0x11, 0x89, 0x00, 0x98  },              //      6
                { 0x11, 0x89, 0x00, 0x98  },              //      7

                { 0x11, 0x90, 0x09, 0x88  },              //      8
                { 0x11, 0x80, 0x09, 0x98  },              //      9
                { 0x11, 0x80, 0x09, 0x98  },              //      10
                { 0x11, 0x90, 0x09, 0x88  },              //      11

                { 0x00, 0x9a, 0xa9, 0x88  },              //      12
                { 0x00, 0x8a, 0xa9, 0x98  },              //      13
                { 0x00, 0x89, 0xaa, 0x98  },              //      14

                { 0x00, 0x9a, 0xa9, 0x88  },              //      15
                { 0x00, 0x8a, 0xa9, 0x98  },              //      16
                { 0x00, 0x89, 0xaa, 0x98  },              //      17

                { 0xbb, 0x89, 0xaa, 0x98  },              //      18
                { 0xbb, 0x8a, 0xa9, 0x98  },              //      19
                { 0xbb, 0x9a, 0xa9, 0x88  },              //      20
                { 0xbb, 0x9a, 0xa9, 0x88  },              //      21
                { 0xbb, 0x8a, 0xa9, 0x98  },              //      22
                { 0xbb, 0x89, 0xaa, 0x98  },              //      23


        // single excitations
                { 0x00, 0x00, 0x00, 0x00  },              //      24
                { 0x00, 0x00, 0x00, 0x00  },              //      25
                { 0x00, 0x00, 0x00, 0x00  },              //      26
                { 0x00, 0x00, 0x00, 0x00  },              //      27


        // no excitation
                { 0x00, 0x00, 0x00, 0x00  }               //      28
	};	
//	addressing:
//	[dK][P][R]
static INT	caseNr[3][5][6] =
	{
		{
			{	6,	8,	9,	10,	11,	7	},
			{	2,	-1,	-1,	-1,	-1,	-1	},
			{	24,	25,	-1,	-1,	-1,	-1	},
			{	1,	-1,	-1,	-1,	-1,	-1	},
			{	28,	-1,	-1,	-1,	-1,	-1	}
		},
		{
			{	15,	16,	17,	12,	13,	14	},
			{	5,	4,	-1,	-1,	-1,	-1	},
			{	26,	27,	-1,	-1,	-1,	-1	},
			{	-1,	-1,	-1,	-1,	-1,	-1	},
			{	-1,	-1,	-1,	-1,	-1,	-1	}
		},
		{
			{	18,	19,	20,	21,	22,	23	},
			{	-1,	-1,	-1,	-1,	-1,	-1	},
			{	-1,	-1,	-1,	-1,	-1,	-1	},
			{	-1,	-1,	-1,	-1,	-1,	-1	},
			{	-1,	-1,	-1,	-1,	-1,	-1	}
		}
	};
		


template <class TMOType>
TableCase<TMOType>::TableCase(
	INT _openShells, INT fall)
{
	openShells = _openShells;
	dK = Cases[fall].dK;
	P = Cases[fall].P;
	R = Cases[fall].R;
	qR = 1;
	qL = 1;
}


template <class TMOType>
const char	*TableCase<TMOType>::getName() const
{
	if ( P==0 || R==0 )
		return "exc>2";
	else
		return Cases[caseNr[abs(dK)][P-1][R-1]].name;
}



template <class TMOType>
void	TableCase<TMOType>::setNr(INT _openShells, INT Nr, INT _qR, INT _qL)
{
	openShells = _openShells;
	dK = Cases[Nr].dK;
	P = Cases[Nr].P;
	R = Cases[Nr].R;
	qR = _qR;
	qL = _qL;
}




template <class TMOType>
INT	TableCase<TMOType>::calc(
	const Configuration<TMOType> & a, const Configuration<TMOType> & b)
{
DiffConf<TMOType>	dc;

	dc.calcDiffConf(a, b);
		
	return calc(dc);
}


template <>
INT	TableCase<MOType>::calc(const DiffConf<MOType> & dc)
{
	dK = (dc.getTo().getNumberOfOpenShells()
		 - dc.getFrom().getNumberOfOpenShells()) >> 1; 

//	cout << "dc = " << dc << endl;



INT	fall = 0;
#define CALC \
	{\
	qL = qR = 0;\
	if ( dc.getFrom().getNumberOfOpenShells() )\
		qL = calcQfromPos(dc.getPosFromP(), dc.getFrom().getNumberOfOpenShells());\
\
	if ( dc.getTo().getNumberOfOpenShells() )\
		qR = calcQfromPos(dc.getPosToP(), dc.getTo().getNumberOfOpenShells());\
\
\
INT	key =	FROM.getNumberOfOpenShells() |\
			(TO.getNumberOfOpenShells() << 2) |\
			(FROM.getNumberOfClosedShells() << 4) |\
			(TO.getNumberOfClosedShells() << 6);\
\
\
	switch ( key ) {\
	case 80:\
		fall = 1;\
		break;\
\
	case 85:\
		if (	FROM.getOpenShell(0)==TO.getClosedShell(0) && \
				TO.getOpenShell(0)==FROM.getClosedShell(0) )\
			fall = 25;\
		else\
		if (	FROM.getOpenShell(0)==TO.getClosedShell(0)  )\
			fall = 2;\
		else\
			fall = 3;\
		break;\
\
	case 24:\
		if ( FROM.getClosedShell(0)==TO.getOpenShell(0) )\
			fall = 26;\
		else\
		if ( FROM.getClosedShell(0)==TO.getOpenShell(1) )\
			fall = 27;\
		else\
			fall = 4;\
		break;\
\
	case 104:\
		if ( 	FROM.getClosedShell(0)==TO.getOpenShell(0) &&\
				FROM.getClosedShell(1)==TO.getOpenShell(1) )\
			fall = 5;\
		break;\
\
	case 10:\
		fall = 6;\
		break;\
\
	case 170:\
		if (	FROM.getClosedShell(0)==TO.getOpenShell(0) &&\
				FROM.getClosedShell(1)==TO.getOpenShell(1) &&\
				TO.getClosedShell(0)==FROM.getOpenShell(0) &&\
				TO.getClosedShell(1)==FROM.getOpenShell(1) )\
			fall = 7;\
		break;\
\
	case 90:\
		if ( FROM.getOpenShell(0)==TO.getClosedShell(0) )\
		{	if ( TO.getOpenShell(0)==FROM.getClosedShell(0) )\
				fall = 8;\
			else\
			if ( TO.getOpenShell(1)==FROM.getClosedShell(0) )\
				fall = 9;\
		}\
		else\
		if ( FROM.getOpenShell(1)==TO.getClosedShell(0) )\
		{	if ( TO.getOpenShell(0)==FROM.getClosedShell(0) )\
				fall = 10;\
			else\
			if ( TO.getOpenShell(1)==FROM.getClosedShell(0) )\
				fall = 11;\
		}\
		break;\
\
	case 29:\
		if ( TO.getOpenShell(0)==FROM.getClosedShell(0) )\
			fall = 12;\
		else\
		if ( TO.getOpenShell(1)==FROM.getClosedShell(0) )\
			fall = 13;\
		else\
		if ( TO.getOpenShell(2)==FROM.getClosedShell(0) )\
			fall = 14;\
		break;\
\
	case 109:\
		if ( FROM.getOpenShell(0)!=TO.getClosedShell(0) )\
			break;\
		if ( TO.getOpenShell(2)==FROM.getClosedShell(1) )\
		{	if ( TO.getOpenShell(1)==FROM.getClosedShell(0) )\
				fall = 15;\
			else\
			if ( TO.getOpenShell(0)==FROM.getClosedShell(0) )\
				fall = 16;\
		}\
		else\
		if ( 	TO.getOpenShell(0)==FROM.getClosedShell(0) &&\
				TO.getOpenShell(1)==FROM.getClosedShell(1) )\
			fall = 17;\
		break;\
\
	case 48:\
		if ( TO.getOpenShell(0)==FROM.getClosedShell(0) )\
		{	if ( TO.getOpenShell(1)==FROM.getClosedShell(1) )\
				fall = 18;\
			else\
			if ( TO.getOpenShell(2)==FROM.getClosedShell(1) )\
				fall = 19;\
			else\
			if ( TO.getOpenShell(3)==FROM.getClosedShell(1) )\
				fall = 20;\
		}\
		else\
		if ( TO.getOpenShell(1)==FROM.getClosedShell(0) )\
		{	if ( TO.getOpenShell(2)==FROM.getClosedShell(1) )\
				fall = 21;\
			else\
			if ( TO.getOpenShell(3)==FROM.getClosedShell(1) )\
				fall = 22;\
		}\
		if ( 	TO.getOpenShell(2)==FROM.getClosedShell(0) &&\
				TO.getOpenShell(3)==FROM.getClosedShell(1) )\
				fall = 23;\
		break;\
\
	case 5:\
		fall = 24;\
		break;\
\
	case 0:\
		fall = 28;\
		break;\
\
	default:\
		fall = 0;\
		break;\
	}\
\
	openShells = FROMSHELLS;\
	P = Cases[fall].P;\
	R = Cases[fall].R;\
\
INT	code;\
TwoElectronIntegralIndex<MOType> Cb;\
TwoElectronIntegralIndex<MOType> Ex;\
	for ( INT i=0 ; i<4 ; i++ )\
	{	code = CoulombExchangeIntegral[fall][i];\
		Cb[i] = (code & 0x08 ?\
				(code & 0x04 ? \
					TO.getClosedShell(code & 3) :\
					TO.getOpenShell(code & 3) ) :\
				(code & 0x04 ? \
					FROM.getClosedShell(code & 3) :\
					FROM.getOpenShell(code & 3) ) );\
		if ( P==1 ) \
		{	code >>= 4;\
			Ex[i] = (code & 0x08 ?\
					(code & 0x04 ? \
						TO.getClosedShell(code & 3) :\
						TO.getOpenShell(code & 3) ) :\
					(code & 0x04 ? \
						FROM.getClosedShell(code & 3) :\
						FROM.getOpenShell(code & 3) ) );\
		}\
	}\
	if ( P==1 )\
		CbExIndex.set(Cb, Ex);\
	else\
		CbExIndex.set(Cb);\
	}


//*********************************************************************
//	The algorithm for table case calculation requires FROM to have
//	less or equal open shells than TO.
//	Therefor the following macros enable the compiler to generate efficient
//	inlined code for the configuration.getXShell() methods.
//	Other solutions would either require to swap the data members
//	of the DiffConf class or would result in less efficient inlined code.
	if ( dK>=0 )
	{
#define FROM dc.getFrom()
#define FROMP dc.getPosFromP()
#define FROMSHELLS dc.getOpenShellsTo()
#define TO dc.getTo()
#define TOP dc.getPosToP()
	CALC
#undef FROM
#undef FROMP
#undef FROMSHELLS
#undef TO
#undef TOP
	}
	else
	{
#define FROM dc.getTo()
#define FROMP dc.getPosToP()
#define FROMSHELLS dc.getOpenShellsFrom()
#define TO dc.getFrom()
#define TOP dc.getPosFromP()
	CALC
#undef FROM
#undef FROMP
#undef FROMSHELLS
#undef TO
#undef TOP
	}
#undef CALC


	return fall!=0;
/*	for ( INT j=0 ; j<2 ; j++ )
	{	cout << "(";
		for ( INT i=0 ; i<4 ; i++ )
		{	cout << c[i][j];
			if ( i==1 )
				cout << "|";
		}
		cout << ")";
	}
*/
}


template <>
INT	TableCase<GeneralizedMO>::calc(const DiffConf<GeneralizedMO> & dc)
{
	dK = (dc.getTo().getNumberOfOpenShells()
		 - dc.getFrom().getNumberOfOpenShells()) >> 1; 

//	cout << "dc = " << dc << endl;



INT	fall = 0;
#define CALC \
	{\
	qL = qR = 0;\
	if ( dc.getFrom().getNumberOfOpenShells() )\
		qL = calcQfromPos(dc.getPosFromP(), dc.getFrom().getNumberOfOpenShells());\
\
	if ( dc.getTo().getNumberOfOpenShells() )\
		qR = calcQfromPos(dc.getPosToP(), dc.getTo().getNumberOfOpenShells());\
\
\
INT	key =	FROM.getNumberOfOpenShells() |\
			(TO.getNumberOfOpenShells() << 2) |\
			(FROM.getNumberOfClosedShells() << 4) |\
			(TO.getNumberOfClosedShells() << 6);\
\
\
	switch ( key ) {\
	case 80:\
		fall = 1;\
		break;\
\
	case 85:\
		if (	FROM.getOpenShell(0)==TO.getClosedShell(0) && \
				TO.getOpenShell(0)==FROM.getClosedShell(0) )\
			fall = 25;\
		else\
		if (	FROM.getOpenShell(0)==TO.getClosedShell(0)  )\
			fall = 2;\
		else\
			fall = 3;\
		break;\
\
	case 24:\
		if ( FROM.getClosedShell(0)==TO.getOpenShell(0) )\
			fall = 26;\
		else\
		if ( FROM.getClosedShell(0)==TO.getOpenShell(1) )\
			fall = 27;\
		else\
			fall = 4;\
		break;\
\
	case 104:\
		if ( 	FROM.getClosedShell(0)==TO.getOpenShell(0) &&\
				FROM.getClosedShell(1)==TO.getOpenShell(1) )\
			fall = 5;\
		break;\
\
	case 10:\
		fall = 6;\
		break;\
\
	case 170:\
		if (	FROM.getClosedShell(0)==TO.getOpenShell(0) &&\
				FROM.getClosedShell(1)==TO.getOpenShell(1) &&\
				TO.getClosedShell(0)==FROM.getOpenShell(0) &&\
				TO.getClosedShell(1)==FROM.getOpenShell(1) )\
			fall = 7;\
		break;\
\
	case 90:\
		if ( FROM.getOpenShell(0)==TO.getClosedShell(0) )\
		{	if ( TO.getOpenShell(0)==FROM.getClosedShell(0) )\
				fall = 8;\
			else\
			if ( TO.getOpenShell(1)==FROM.getClosedShell(0) )\
				fall = 9;\
		}\
		else\
		if ( FROM.getOpenShell(1)==TO.getClosedShell(0) )\
		{	if ( TO.getOpenShell(0)==FROM.getClosedShell(0) )\
				fall = 10;\
			else\
			if ( TO.getOpenShell(1)==FROM.getClosedShell(0) )\
				fall = 11;\
		}\
		break;\
\
	case 29:\
		if ( TO.getOpenShell(0)==FROM.getClosedShell(0) )\
			fall = 12;\
		else\
		if ( TO.getOpenShell(1)==FROM.getClosedShell(0) )\
			fall = 13;\
		else\
		if ( TO.getOpenShell(2)==FROM.getClosedShell(0) )\
			fall = 14;\
		break;\
\
	case 109:\
		if ( FROM.getOpenShell(0)!=TO.getClosedShell(0) )\
			break;\
		if ( TO.getOpenShell(2)==FROM.getClosedShell(1) )\
		{	if ( TO.getOpenShell(1)==FROM.getClosedShell(0) )\
				fall = 15;\
			else\
			if ( TO.getOpenShell(0)==FROM.getClosedShell(0) )\
				fall = 16;\
		}\
		else\
		if ( 	TO.getOpenShell(0)==FROM.getClosedShell(0) &&\
				TO.getOpenShell(1)==FROM.getClosedShell(1) )\
			fall = 17;\
		break;\
\
	case 48:\
		if ( TO.getOpenShell(0)==FROM.getClosedShell(0) )\
		{	if ( TO.getOpenShell(1)==FROM.getClosedShell(1) )\
				fall = 18;\
			else\
			if ( TO.getOpenShell(2)==FROM.getClosedShell(1) )\
				fall = 19;\
			else\
			if ( TO.getOpenShell(3)==FROM.getClosedShell(1) )\
				fall = 20;\
		}\
		else\
		if ( TO.getOpenShell(1)==FROM.getClosedShell(0) )\
		{	if ( TO.getOpenShell(2)==FROM.getClosedShell(1) )\
				fall = 21;\
			else\
			if ( TO.getOpenShell(3)==FROM.getClosedShell(1) )\
				fall = 22;\
		}\
		if ( 	TO.getOpenShell(2)==FROM.getClosedShell(0) &&\
				TO.getOpenShell(3)==FROM.getClosedShell(1) )\
				fall = 23;\
		break;\
\
	case 5:\
		fall = 24;\
		break;\
\
	case 0:\
		fall = 28;\
		break;\
\
	default:\
		fall = 0;\
		break;\
	}\
\
	openShells = FROMSHELLS;\
	P = Cases[fall].P;\
	R = Cases[fall].R;\
	}


//*********************************************************************
//	The algorithm for table case calculation requires FROM to have
//	less or equal open shells than TO.
//	Therefor the following macros enable the compiler to generate efficient
//	inlined code for the configuration.getXShell() methods.
//	Other solutions would either require to swap the data members
//	of the DiffConf class or would result in less efficient inlined code.
	if ( dK>=0 )
	{
#define FROM dc.getFrom()
#define FROMP dc.getPosFromP()
#define FROMSHELLS dc.getOpenShellsTo()
#define TO dc.getTo()
#define TOP dc.getPosToP()
	CALC
#undef FROM
#undef FROMP
#undef FROMSHELLS
#undef TO
#undef TOP
	}
	else
	{
#define FROM dc.getTo()
#define FROMP dc.getPosToP()
#define FROMSHELLS dc.getOpenShellsFrom()
#define TO dc.getFrom()
#define TOP dc.getPosFromP()
	CALC
#undef FROM
#undef FROMP
#undef FROMSHELLS
#undef TO
#undef TOP
	}
#undef CALC


	return fall!=0;
/*	for ( INT j=0 ; j<2 ; j++ )
	{	cout << "(";
		for ( INT i=0 ; i<4 ; i++ )
		{	cout << c[i][j];
			if ( i==1 )
				cout << "|";
		}
		cout << ")";
	}
*/
}



template <class TMOType>
void	TableCase<TMOType>::calcLess3NoInd(
	const Configuration<TMOType> & a, const Configuration<TMOType> & b)
{
DiffConf<TMOType>	dc;

	dc.calcDiffConf(a, b);
	calcLess3NoInd(dc);
}


template <class TMOType>
void	TableCase<TMOType>::calcLess3NoInd(const DiffConf<TMOType> & dc)
{
	dK = (dc.getTo().getNumberOfOpenShells()
		 - dc.getFrom().getNumberOfOpenShells()) >> 1; 

//	cout << "dc = " << dc << endl;



#define CALC \
	{\
	qL = qR = 0;\
	if ( dc.getFrom().getNumberOfOpenShells() )\
		qL = calcQfromPos(dc.getPosFromP(), dc.getFrom().getNumberOfOpenShells());\
\
	if ( dc.getTo().getNumberOfOpenShells() )\
		qR = calcQfromPos(dc.getPosToP(), dc.getTo().getNumberOfOpenShells());\
\
\
INT	key =	FROM.getNumberOfOpenShells() |\
			(TO.getNumberOfOpenShells() << 2) |\
			(FROM.getNumberOfClosedShells() << 4) |\
			(TO.getNumberOfClosedShells() << 6);\
\
INT	fall;\
	switch ( key ) {\
	case 80:\
		fall = 1;\
		break;\
\
	case 85:\
		if (	FROM.getOpenShell(0)==TO.getClosedShell(0) && \
				TO.getOpenShell(0)==FROM.getClosedShell(0) )\
			fall = 25;\
		else\
		if (	FROM.getOpenShell(0)==TO.getClosedShell(0) )\
			fall = 2;\
		else\
			fall = 3;\
		break;\
\
	case 24:\
		if ( FROM.getClosedShell(0)==TO.getOpenShell(0) )\
			fall = 26;\
		else\
		if ( FROM.getClosedShell(0)==TO.getOpenShell(1) )\
			fall = 27;\
		else\
			fall = 4;\
		break;\
\
	case 104:\
		fall = 5;\
		break;\
\
	case 10:\
		fall = 6;\
		break;\
\
	case 170:\
		fall = 7;\
		break;\
\
	case 90:\
		if ( FROM.getOpenShell(0)==TO.getClosedShell(0) )\
		{	if ( TO.getOpenShell(0)==FROM.getClosedShell(0) )\
				fall = 8;\
			else\
				fall = 9;\
		}\
		else\
		{	if ( TO.getOpenShell(0)==FROM.getClosedShell(0) )\
				fall = 10;\
			else\
				fall = 11;\
		}\
		break;\
\
	case 29:\
		if ( TO.getOpenShell(0)==FROM.getClosedShell(0) )\
			fall = 12;\
		else\
		if ( TO.getOpenShell(1)==FROM.getClosedShell(0) )\
			fall = 13;\
		else\
			fall = 14;\
		break;\
\
	case 109:\
		if ( TO.getOpenShell(1)==FROM.getClosedShell(0) )\
			fall = 15;\
		else\
		if ( TO.getOpenShell(1)==FROM.getClosedShell(1) )\
			fall = 17;\
		else\
			fall = 16;\
		break;\
\
	case 48:\
		if ( TO.getOpenShell(0)==FROM.getClosedShell(0) )\
		{	if ( TO.getOpenShell(1)==FROM.getClosedShell(1) )\
				fall = 18;\
			else\
			if ( TO.getOpenShell(2)==FROM.getClosedShell(1) )\
				fall = 19;\
			else\
				fall = 20;\
		}\
		else\
		{	if ( TO.getOpenShell(2)==FROM.getClosedShell(1) )\
				fall = 21;\
			else\
			if ( TO.getOpenShell(2)==FROM.getClosedShell(0) )\
				fall = 23;\
			else\
				fall = 22;\
		}\
		break;\
\
	case 5:\
		fall = 24;\
		break;\
\
	case 0:\
		fall = 28;\
		break;\
\
	default:\
		fall = 0;\
		break;\
	}\
\
	openShells = FROMSHELLS;\
	P = Cases[fall].P;\
	R = Cases[fall].R;\
	}

//*********************************************************************
//	The algorithm for table case calculation requires FROM to have
//	less or equal open shells than TO.
//	Therefor the following macros enable the compiler to generate efficient
//	inlined code for the configuration.getXShell() methods.
//	Other solutions would either require to swap the data members
//	of the DiffConf class or would result in less efficient inlined code.
	if ( dK>=0 )
	{
#define FROM dc.getFrom()
#define FROMP dc.getPosFromP()
#define FROMSHELLS dc.getOpenShellsTo()
#define TO dc.getTo()
#define TOP dc.getPosToP()
	CALC
#undef FROM
#undef FROMP
#undef FROMSHELLS
#undef TO
#undef TOP
	}
	else
	{
#define FROM dc.getTo()
#define FROMP dc.getPosToP()
#define FROMSHELLS dc.getOpenShellsFrom()
#define TO dc.getFrom()
#define TOP dc.getPosFromP()
	CALC
#undef FROM
#undef FROMP
#undef FROMSHELLS
#undef TO
#undef TOP
	}
#undef CALC

//	cout << *this << endl;

}



template <>
void	TableCase<MOType>::calcLess3(const DiffConf<MOType> & dc)
{
	dK = (dc.getTo().getNumberOfOpenShells()
		 - dc.getFrom().getNumberOfOpenShells()) >> 1; 
//	cout << "dc = " << dc << endl;



#define CALC \
	{\
	qL = qR = 0;\
	if ( dc.getFrom().getNumberOfOpenShells() )\
		qL = calcQfromPos(dc.getPosFromP(), dc.getFrom().getNumberOfOpenShells());\
\
	if ( dc.getTo().getNumberOfOpenShells() )\
		qR = calcQfromPos(dc.getPosToP(), dc.getTo().getNumberOfOpenShells());\
\
\
INT	key =	FROM.getNumberOfOpenShells() |\
			(TO.getNumberOfOpenShells() << 2) |\
			(FROM.getNumberOfClosedShells() << 4) |\
			(TO.getNumberOfClosedShells() << 6);\
\
INT	fall;\
	switch ( key ) {\
	case 80:\
		fall = 1;\
		break;\
\
	case 85:\
		if (	FROM.getOpenShell(0)==TO.getClosedShell(0) && \
				TO.getOpenShell(0)==FROM.getClosedShell(0) )\
			fall = 25;\
		else\
		if (	FROM.getOpenShell(0)==TO.getClosedShell(0) )\
			fall = 2;\
		else\
			fall = 3;\
		break;\
\
	case 24:\
		if ( FROM.getClosedShell(0)==TO.getOpenShell(0) )\
			fall = 26;\
		else\
		if ( FROM.getClosedShell(0)==TO.getOpenShell(1) )\
			fall = 27;\
		else\
			fall = 4;\
		break;\
\
	case 104:\
		fall = 5;\
		break;\
\
	case 10:\
		fall = 6;\
		break;\
\
	case 170:\
		fall = 7;\
		break;\
\
	case 90:\
		if ( FROM.getOpenShell(0)==TO.getClosedShell(0) )\
		{	if ( TO.getOpenShell(0)==FROM.getClosedShell(0) )\
				fall = 8;\
			else\
				fall = 9;\
		}\
		else\
		{	if ( TO.getOpenShell(0)==FROM.getClosedShell(0) )\
				fall = 10;\
			else\
				fall = 11;\
		}\
		break;\
\
	case 29:\
		if ( TO.getOpenShell(0)==FROM.getClosedShell(0) )\
			fall = 12;\
		else\
		if ( TO.getOpenShell(1)==FROM.getClosedShell(0) )\
			fall = 13;\
		else\
			fall = 14;\
		break;\
\
	case 109:\
		if ( TO.getOpenShell(1)==FROM.getClosedShell(0) )\
			fall = 15;\
		else\
		if ( TO.getOpenShell(1)==FROM.getClosedShell(1) )\
			fall = 17;\
		else\
			fall = 16;\
		break;\
\
	case 48:\
		if ( TO.getOpenShell(0)==FROM.getClosedShell(0) )\
		{	if ( TO.getOpenShell(1)==FROM.getClosedShell(1) )\
				fall = 18;\
			else\
			if ( TO.getOpenShell(2)==FROM.getClosedShell(1) )\
				fall = 19;\
			else\
				fall = 20;\
		}\
		else\
		{	if ( TO.getOpenShell(2)==FROM.getClosedShell(1) )\
				fall = 21;\
			else\
			if ( TO.getOpenShell(2)==FROM.getClosedShell(0) )\
				fall = 23;\
			else\
				fall = 22;\
		}\
		break;\
\
	case 5:\
		fall = 24;\
		break;\
\
	case 0:\
		fall = 28;\
		break;\
\
	default:\
		fall = 0;\
		break;\
	}\
\
	openShells = FROMSHELLS;\
	P = Cases[fall].P;\
	R = Cases[fall].R;\
\
	if ( P==3 || P==5 )\
		return;\
\
INT	code;\
TwoElectronIntegralIndex<MOType> Cb;\
TwoElectronIntegralIndex<MOType> Ex;\
	for ( INT i=0 ; i<4 ; i++ )\
	{	code = CoulombExchangeIntegral[fall][i];\
		Cb[i] = (code & 0x08 ?\
				(code & 0x04 ? \
					TO.getClosedShell(code & 3) :\
					TO.getOpenShell(code & 3) ) :\
				(code & 0x04 ? \
					FROM.getClosedShell(code & 3) :\
					FROM.getOpenShell(code & 3) ) );\
		if ( P==1 ) \
		{	code >>= 4;\
			Ex[i] = (code & 0x08 ?\
					(code & 0x04 ? \
						TO.getClosedShell(code & 3) :\
						TO.getOpenShell(code & 3) ) :\
					(code & 0x04 ? \
						FROM.getClosedShell(code & 3) :\
						FROM.getOpenShell(code & 3) ) );\
		}\
	}\
	if ( P==1 )\
		CbExIndex.set(Cb, Ex);\
	else\
		CbExIndex.set(Cb);\
	}


//*********************************************************************
//	The algorithm for table case calculation requires FROM to have
//	less or equal open shells than TO.
//	Therefor the following macros enable the compiler to generate efficient
//	inlined code for the configuration.getXShell() methods.
//	Other solutions would either require to swap the data members
//	of the DiffConf class or would result in less efficient inlined code.
	if ( dK>=0 )
	{
#define FROM dc.getFrom()
#define FROMP dc.getPosFromP()
#define FROMSHELLS dc.getOpenShellsTo()
#define TO dc.getTo()
#define TOP dc.getPosToP()
	CALC
#undef FROM
#undef FROMP
#undef FROMSHELLS
#undef TO
#undef TOP
	}
	else
	{
#define FROM dc.getTo()
#define FROMP dc.getPosToP()
#define FROMSHELLS dc.getOpenShellsFrom()
#define TO dc.getFrom()
#define TOP dc.getPosFromP()
	CALC
#undef FROM
#undef FROMP
#undef FROMSHELLS
#undef TO
#undef TOP
	}
#undef CALC

//	cout << *this << endl;
//	qR = 6;

}


template <class TMOType>
void	TableCase<TMOType>::calcLess3(const DiffConf<TMOType> & dc,
		TwoElectronIntegralIndex<TMOType> & Cb,
		TwoElectronIntegralIndex<TMOType> & Ex)
{
	dK = (dc.getTo().getNumberOfOpenShells()
		 - dc.getFrom().getNumberOfOpenShells()) >> 1; 

//	cout << "dc = " << dc << endl;



#define CALC \
	{\
	qL = qR = 0;\
	if ( dc.getFrom().getNumberOfOpenShells() )\
		qL = calcQfromPos(dc.getPosFromP(), dc.getFrom().getNumberOfOpenShells());\
\
	if ( dc.getTo().getNumberOfOpenShells() )\
		qR = calcQfromPos(dc.getPosToP(), dc.getTo().getNumberOfOpenShells());\
\
\
INT	key =	FROM.getNumberOfOpenShells() |\
			(TO.getNumberOfOpenShells() << 2) |\
			(FROM.getNumberOfClosedShells() << 4) |\
			(TO.getNumberOfClosedShells() << 6);\
\
INT	fall;\
	switch ( key ) {\
	case 80:\
		fall = 1;\
		break;\
\
	case 85:\
		if (	FROM.getOpenShell(0)==TO.getClosedShell(0) && \
				TO.getOpenShell(0)==FROM.getClosedShell(0) )\
			fall = 25;\
		else\
		if (	FROM.getOpenShell(0)==TO.getClosedShell(0) )\
			fall = 2;\
		else\
			fall = 3;\
		break;\
\
	case 24:\
		if ( FROM.getClosedShell(0)==TO.getOpenShell(0) )\
			fall = 26;\
		else\
		if ( FROM.getClosedShell(0)==TO.getOpenShell(1) )\
			fall = 27;\
		else\
			fall = 4;\
		break;\
\
	case 104:\
		fall = 5;\
		break;\
\
	case 10:\
		fall = 6;\
		break;\
\
	case 170:\
		fall = 7;\
		break;\
\
	case 90:\
		if ( FROM.getOpenShell(0)==TO.getClosedShell(0) )\
		{	if ( TO.getOpenShell(0)==FROM.getClosedShell(0) )\
				fall = 8;\
			else\
				fall = 9;\
		}\
		else\
		{	if ( TO.getOpenShell(0)==FROM.getClosedShell(0) )\
				fall = 10;\
			else\
				fall = 11;\
		}\
		break;\
\
	case 29:\
		if ( TO.getOpenShell(0)==FROM.getClosedShell(0) )\
			fall = 12;\
		else\
		if ( TO.getOpenShell(1)==FROM.getClosedShell(0) )\
			fall = 13;\
		else\
			fall = 14;\
		break;\
\
	case 109:\
		if ( TO.getOpenShell(1)==FROM.getClosedShell(0) )\
			fall = 15;\
		else\
		if ( TO.getOpenShell(1)==FROM.getClosedShell(1) )\
			fall = 17;\
		else\
			fall = 16;\
		break;\
\
	case 48:\
		if ( TO.getOpenShell(0)==FROM.getClosedShell(0) )\
		{	if ( TO.getOpenShell(1)==FROM.getClosedShell(1) )\
				fall = 18;\
			else\
			if ( TO.getOpenShell(2)==FROM.getClosedShell(1) )\
				fall = 19;\
			else\
				fall = 20;\
		}\
		else\
		{	if ( TO.getOpenShell(2)==FROM.getClosedShell(1) )\
				fall = 21;\
			else\
			if ( TO.getOpenShell(2)==FROM.getClosedShell(0) )\
				fall = 23;\
			else\
				fall = 22;\
		}\
		break;\
\
	case 5:\
		fall = 24;\
		break;\
\
	case 0:\
		fall = 28;\
		break;\
\
	default:\
		fall = 0;\
		break;\
	}\
\
	openShells = FROMSHELLS;\
	P = Cases[fall].P;\
	R = Cases[fall].R;\
\
	if ( P==3 || P==5 )\
		return;\
\
INT	code;\
	for ( INT i=0 ; i<4 ; i++ )\
	{	code = CoulombExchangeIntegral[fall][i];\
		Cb[i] = (code & 0x08 ?\
				(code & 0x04 ? \
					TO.getClosedShell(code & 3) :\
					TO.getOpenShell(code & 3) ) :\
				(code & 0x04 ? \
					FROM.getClosedShell(code & 3) :\
					FROM.getOpenShell(code & 3) ) );\
		if ( P==1 ) \
		{	code >>= 4;\
			Ex[i] = (code & 0x08 ?\
					(code & 0x04 ? \
						TO.getClosedShell(code & 3) :\
						TO.getOpenShell(code & 3) ) :\
					(code & 0x04 ? \
						FROM.getClosedShell(code & 3) :\
						FROM.getOpenShell(code & 3) ) );\
		}\
	}\
	}


//*********************************************************************
//	The algorithm for table case calculation requires FROM to have
//	less or equal open shells than TO.
//	Therefor the following macros enable the compiler to generate efficient
//	inlined code for the configuration.getXShell() methods.
//	Other solutions would either require to swap the data members
//	of the DiffConf class or would result in less efficient inlined code.
	if ( dK>=0 )
	{
#define FROM dc.getFrom()
#define FROMP dc.getPosFromP()
#define FROMSHELLS dc.getOpenShellsTo()
#define TO dc.getTo()
#define TOP dc.getPosToP()
	CALC
#undef FROM
#undef FROMP
#undef FROMSHELLS
#undef TO
#undef TOP
	}
	else
	{
#define FROM dc.getTo()
#define FROMP dc.getPosToP()
#define FROMSHELLS dc.getOpenShellsFrom()
#define TO dc.getFrom()
#define TOP dc.getPosFromP()
	CALC
#undef FROM
#undef FROMP
#undef FROMSHELLS
#undef TO
#undef TOP
	}
#undef CALC

//	cout << *this << endl;
//	qR = 6;

}



template <class TMOType>
void	TableCase<TMOType>::calcQ(
	const Configuration<TMOType> & a, const Configuration<TMOType> & b)
{
	dK = (b.getNumberOfOpenShells() - a.getNumberOfOpenShells()) >> 1; 

DiffConf<TMOType>	dc;

	if ( dK>=0 )
		dc.calcDiffConf(a, b);
	else
		dc.calcDiffConf(b, a);

	qL = qR = 0;
	if ( dc.getFrom().getNumberOfOpenShells() )
		qL = calcQfromPos(dc.getPosFromP(), dc.getFrom().getNumberOfOpenShells());

	if ( dc.getTo().getNumberOfOpenShells() )
		qR = calcQfromPos(dc.getPosToP(), dc.getTo().getNumberOfOpenShells());
}


template <class TMOType>
void	TableCase<TMOType>::calcQ(const DiffConf<TMOType> & dc)
{
//	cout << "dc = " << dc << endl;
	dK = (dc.getTo().getNumberOfOpenShells()
		 - dc.getFrom().getNumberOfOpenShells()) >> 1; 

	qL = qR = 0;
	if ( dc.getFrom().getNumberOfOpenShells() )
		qL = calcQfromPos(dc.getPosFromP(), dc.getFrom().getNumberOfOpenShells());

	if ( dc.getTo().getNumberOfOpenShells() )
		qR = calcQfromPos(dc.getPosToP(), dc.getTo().getNumberOfOpenShells());
}


		

//-----------------------------------------------------------------------

template <class TMOType>
ostream& operator<<(ostream& s, const TableCase<TMOType> & tc)
{
	switch ( tc.getOutputMode() ) {
	case MathObject::TerminalASCII:
/*		s << tc.getOpenShells();
		s << ", " << tc.getName();
		if ( tc.getqR()>0 )
			s << " " << tc.getqR();
		if ( tc.getqL()>0 )
			s << " " << tc.getqL();
		break;
*/
		
		s << "open shells= " << tc.getNumberOfMoreOpenShells();
		s << ", Fall Nr. " << tc.getName() << ": ";
		s << "deltaK=" << tc.getdK();
		s << ", P=" << tc.getP();
		s << ", R=" << tc.getR();
		if ( tc.getqL()>0 )
			s << ", qL=" << tc.getqL();
		if ( tc.getqR()>0 )
			s << ", qR=" << tc.getqR();
		s << tc.CbExIndex << endl;
		break;
		
	case MathObject::TeX:
		break;
	}
	return s;
}

template class TableCase<MOType>;
template class TableCase<GeneralizedMO>;

template ostream& operator << (ostream& s, const
    TableCase<GeneralizedMO> & v);
template ostream& operator << (ostream& s, const
    TableCase<MOType> & v);


//**********************************************************************
//	explicit instantiation for g++ 2.7.2 
//	("template <class TMOType> void calcInteraction(...)" does not work)
/*
void	TableCaseInstanciateTemplates()
{
#define TEMPLATE_INSTANCIATE \
	{\
	BinomialCoefficient	*bin = NULL;\
	Configuration<T> a, b;\
	TableCase<T> t(bin);\
		t.calc(a, b);\
		cout << t;\
	}

#define T MOType
TEMPLATE_INSTANCIATE
#undef T	
	
#define T GeneralizedMO
TEMPLATE_INSTANCIATE
#undef T	
	
#undef TEMPLATE_INSTANCIATE
}
*/

