//***********************************************************************
//
//	Name:			MOIrReps.cc
//
//	Description:	
//
//	Author:			Michael Hanrath
//
//	Version:		0.0
//
//	Date:			12.08.1996
//
//
//
//
//
//***********************************************************************

#include <iostream>

#include "MOIrReps.h"

#include "../IO/Fortran/Fort31FirstRecord.h"

using namespace std;

const INT MOIrReps::maxIrreps = 8;
const IrRep MOIrReps::prodTab[] = {
	0, 1, 2, 3, 4, 5, 6, 7,
	1, 0, 3, 2, 5, 4, 7, 6,
	2, 3, 0, 1, 6, 7, 4, 5,
	3, 2, 1, 0, 7, 6, 5, 4,
	4, 5, 6, 7, 0, 1, 2, 3,
	5, 4, 7, 6, 1, 0, 3, 2,
	6, 7, 4, 5, 2, 3, 0, 1,
	7, 6, 5, 4, 3, 2, 1, 0
	};



MOIrReps::MOIrReps()
{
	maxMO = 0;
	IrReps = 0;
	inIrRep = NULL;
	MOSymmetry = NULL;
	MOinIrRep = NULL;
	IrRepStart = NULL;
//	prodTab = NULL;
}

MOIrReps::MOIrReps(istream &is)
{
	is >> IrReps;
	is >> maxMO;

	inIrRep = new INT [IrReps];
	for ( INT i=0 ; i<IrReps ; i++ )
		is >> inIrRep[i];

/*	prodTab = new IrRep[maxIrreps*maxIrreps];
	for ( INT i=0 ; i<IrReps ; i++ )
		for ( INT j=0 ; j<IrReps ; j++ )
			is >> getProd(i, j);
*/
	init();
}



MOIrReps::MOIrReps(Fort31File _f31)
{
	

	switch ( _f31.format ) {
	case RIFormat:
		{


		//	ACHTUNG:
		//	Produkttafel gilt nur, wenn die Irreps in turbomol
		//	so geordnet sind, dass 
		//  ir1 und ir2 eine Untergruppe bilden
		//	ir1 bis ir4 die naechst groessere
		//	ir1 bis ir8 dann d2h bilden
		//	eventuell Erweiterung noetig, falls
		//	turbomol von dieser Konvention abweicht.
		//	? vielleicht mit map Vertauschung entsprechender
		//	Zeilen und Reihen vornehmen
		//	Standardreihenfolge fuer jede Symmetriegruppe fest vorgeben
		//	eingelesene Reihenfolge damit vergleichen
		//	Reihen und Zeilen entsprechend vertauschen
		//	weitererer Nachteil: es muessten weitere Files eingelesen werden



  			RISym f31(string(_f31.name));
  			IrReps = f31.getN();
  			inIrRep = new INT [IrReps];
  			for ( INT i=0 ; i<IrReps ; i++ )
    				inIrRep[i] = f31.getLJ(i);

		/*  	prodTab = new IrRep[maxIrreps*maxIrreps];
  			INT	k = 0;
  			for ( INT i=0 ; i<IrReps ; i++ )
    				for ( INT j=0 ; j<=i ; j++ )
      					//      prodTab[j*IrReps + i] = prodTab[i*IrReps + j] = f31.getJAB(k++) - 1;
      					prodTab[j*maxIrreps + i] = prodTab[i*maxIrreps + j] = f31.getJAB(i*IrReps +j) - 1;

		*/	init();
	}
		break;

	case GAMESSC1:
		{
		INT	n;
		FortranFileIO	fort31(_f31.name);
			fort31.read(&n);
			maxMO = n;
			IrReps = 1;
			inIrRep = new INT [IrReps];
			inIrRep[0] = maxMO;
			init();
		}
		break;
	
	default:
		{
			Fort31FirstRecord	f31(_f31);


			maxMO = 0;
			IrReps = f31.getN();
			inIrRep = new INT [IrReps];
			for ( INT i=0 ; i<IrReps ; i++ )
				inIrRep[i] = f31.getLJ(i);


		/*	prodTab = new IrRep[maxIrreps*maxIrreps];
		INT	k = 0;



			for ( INT i=0 ; i<maxIrreps ; i++ )
				for ( INT j=0 ; j<=i ; j++ )
					prodTab[j*maxIrreps + i] = prodTab[i*maxIrreps + j] = f31.getJAB(k++) - 1;
		*/
			init();
		}
	}
}

void	MOIrReps::init()
{
	maxMO = 0;
	for ( INT i=0 ; i<IrReps ; i++ )
		maxMO += inIrRep[i];


	MOSymmetry = new INT [maxMO+1];
INT	k = 1;
	for ( INT i=0 ; i<IrReps ; i++ )
		for ( INT j=0 ; j<inIrRep[i] && k<=maxMO ; j++ )
			MOSymmetry[k++] = i;


	MOinIrRep = new MOType[maxMO+1];
	IrRepStart = new MOType[IrReps];

INT sum = 0;
	IrRepStart[0] = 1;
	for ( INT i=0 ; i<IrReps ; i++ )
	{	for ( INT j=0 ; j<inIrRep[i] ; j++ )
			MOinIrRep[j+sum] = j;

		if ( i )
			IrRepStart[i] = sum + 1;
		sum += inIrRep[i];
	}
}


MOIrReps::~MOIrReps()
{
	if ( inIrRep )
		delete[] inIrRep;
	if ( MOSymmetry )
		delete[] MOSymmetry;
	if ( MOinIrRep )
		delete[] MOinIrRep;
	if ( IrRepStart )
		delete[] IrRepStart;
//	if ( prodTab )
//		delete[] prodTab;
}


MOIrReps::MOIrReps(const MOIrReps & irreps)
{
	maxMO = irreps.getMaxMO();
	IrReps = irreps.getNumberOfIrReps();

	inIrRep = new INT [IrReps];
	memcpy(inIrRep, irreps.getInIrRepP(), IrReps*sizeof(INT));

	MOSymmetry = new INT [maxMO+1];
	memcpy(MOSymmetry, irreps.getMOSymmetryP(), (maxMO+1)*sizeof(IrRep));

	MOinIrRep = new MOType[maxMO+1];
	memcpy(MOinIrRep, irreps.getNumberInIrRepP(), (maxMO+1)*sizeof(MOType));

	IrRepStart = new MOType[IrReps];
	memcpy(IrRepStart, irreps.getStartMOP(), IrReps*sizeof(MOType));

/*	if ( irreps.getProdTabP() )
	{
		prodTab = new IrRep[maxIrreps*maxIrreps];
		memcpy(prodTab, irreps.getProdTabP(), maxIrreps*maxIrreps*sizeof(IrRep));
	}
	else
		prodTab = NULL;
*/
}


MOIrReps & MOIrReps::operator = (const MOIrReps & irreps)
{
	maxMO = irreps.getMaxMO();
	IrReps = irreps.getNumberOfIrReps();

	if ( inIrRep )
		delete[] inIrRep;
	inIrRep = new INT [IrReps];
	memcpy(inIrRep, irreps.getInIrRepP(), IrReps*sizeof(INT));

	if ( MOSymmetry )
		delete[] MOSymmetry;
	MOSymmetry = new INT [maxMO+1];
	memcpy(MOSymmetry, irreps.getMOSymmetryP(), (maxMO+1)*sizeof(IrRep));

	if ( MOinIrRep )
		delete[] MOinIrRep;
	MOinIrRep = new MOType[maxMO+1];
	memcpy(MOinIrRep, irreps.getNumberInIrRepP(), (maxMO+1)*sizeof(MOType));

	if ( IrRepStart )
		delete[] IrRepStart;
	IrRepStart = new MOType[IrReps];
	memcpy(IrRepStart, irreps.getStartMOP(), IrReps*sizeof(MOType));

/*	if ( prodTab )
		delete[] prodTab;
	
	if ( irreps.getProdTabP() )
	{
		prodTab = new IrRep[maxIrreps*maxIrreps];
		memcpy(prodTab, irreps.getProdTabP(), maxIrreps*maxIrreps*sizeof(IrRep));
	}
	else
		prodTab = NULL;
*/
	return *this;
}

ostream & MOIrReps::writeToStream(ostream &os)
{
	os << IrReps << endl;
	os << getMaxMO() << endl;
	for ( INT i=0 ; i<getNumberOfIrReps() ; i++ )
		os << getInIrRep(i) << " ";
	os << endl;
/*	for ( INT i=0 ; i<getNumberOfIrReps() ; i++ )
	{	for ( INT j=0 ; j<getNumberOfIrReps() ; j++ )
//			os << (prodTab ? getProd(i, j) : 0 ) << " ";
			os << getProd(i, j)  << " ";
		os << endl;
	}
*/
	return os;	
}


ostream& operator<<(ostream & s, const MOIrReps & moirreps)
{
	s << "maxMO = " << moirreps.getMaxMO() << endl;
	s << endl;
	s << "MOs per symmetry:" << endl;
	for ( INT i=0 ; i<moirreps.getNumberOfIrReps() ; i++ )
		s << moirreps.getInIrRep(i) << " ";
	s << endl;
	s << endl;
	s << "product table:" << endl;
	for ( INT i=0 ; i<8+0*moirreps.getNumberOfIrReps() ; i++ )
	{	for ( INT j=0 ; j<8+0*moirreps.getNumberOfIrReps() ; j++ )
			s << moirreps.getProd(i, j) << " ";
		s << endl;
	}
	s << endl;
	
	return s;
}
