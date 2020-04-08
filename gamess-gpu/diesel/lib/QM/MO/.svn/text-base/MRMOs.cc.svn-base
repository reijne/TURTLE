//***********************************************************************
//
//	Name:			MRMOs.cc
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

#include "MRMOs.h"

#include <string>
#include "MOMapping.h"

using namespace std;

MRMOs::MRMOs() : 
	MOIrReps()
{
	isIntern = NULL;
	inIrRepInt = NULL;
	inIrRepExt = NULL;
	MOinIrRepIntExt = NULL;
	MOIrRepIntExtStart = NULL;
	nInactiveMOs = 0;
	inactiveMO = NULL;
}

MRMOs::MRMOs(Fort31File f31) : 
	MOIrReps(f31)
{
	isIntern = new unsigned char [getMaxMO()+1];
	memset(isIntern, 0, (getMaxMO()+1)*sizeof(unsigned char));
	inIrRepInt = NULL;
	inIrRepExt = NULL;
	MOinIrRepIntExt = NULL;
	MOIrRepIntExtStart = NULL;
	nInactiveMOs = 0;
	inactiveMO = NULL;
}


MRMOs::MRMOs(const MRMOs & mrmos) :
	MOIrReps(mrmos)
{
	isIntern = new unsigned char [getMaxMO()+1];
	memcpy(isIntern, mrmos.getInternP(), (getMaxMO()+1)*sizeof(unsigned char));

	inIrRepInt = new INT [IrReps];
	if ( mrmos.getInIrRepIntP() )
		memcpy(inIrRepInt, mrmos.getInIrRepIntP(), IrReps*sizeof(INT));

	inIrRepExt = new INT [IrReps];
	if ( mrmos.getInIrRepExtP() )
		memcpy(inIrRepExt, mrmos.getInIrRepExtP(), IrReps*sizeof(INT));

	MOinIrRepIntExt = new MOType[maxMO+1];
	if ( mrmos.getNumberInIrRepIntExtP() )
		memcpy(MOinIrRepIntExt, mrmos.getNumberInIrRepIntExtP(),
			(maxMO+1)*sizeof(MOType));

	MOIrRepIntExtStart = new INT [IrReps*2];
	if ( mrmos.getIrRepIntExtStartP() )
		memcpy(MOIrRepIntExtStart, mrmos.getIrRepIntExtStartP(),
			2*IrReps*sizeof(INT));

	nInactiveMOs = mrmos.getNumberOfInactiveMOs();
	inactiveMO = new MOType[nInactiveMOs];
	memcpy(inactiveMO, mrmos.getInactiveMOP(), nInactiveMOs*sizeof(MOType));

}


MRMOs & MRMOs::operator = (const MRMOs & mrmos)
{
//	*((MOIrReps *) this) = *((MOIrReps *) &mrmos);
	(MOIrReps &) *this = mrmos;

	if ( isIntern )
		delete[] isIntern;
	isIntern = new unsigned char [getMaxMO()+1];
	memcpy(isIntern, mrmos.getInternP(), (getMaxMO()+1)*sizeof(unsigned char));

	if ( inIrRepInt )
		delete[] inIrRepInt;
	inIrRepInt = new INT [IrReps];
	if ( mrmos.getInIrRepIntP() )
		memcpy(inIrRepInt, mrmos.getInIrRepIntP(), IrReps*sizeof(INT));

	if ( inIrRepExt )
		delete[] inIrRepExt;
	inIrRepExt = new INT [IrReps];
	if ( mrmos.getInIrRepExtP() )
		memcpy(inIrRepExt, mrmos.getInIrRepExtP(), IrReps*sizeof(INT));

	if ( MOinIrRepIntExt )
		delete[] MOinIrRepIntExt;
	MOinIrRepIntExt = new MOType[maxMO+1];
	if ( mrmos.getNumberInIrRepIntExtP() )
		memcpy(MOinIrRepIntExt, mrmos.getNumberInIrRepIntExtP(),
			(maxMO+1)*sizeof(MOType));

	if ( MOIrRepIntExtStart )
		delete[] MOIrRepIntExtStart;
	MOIrRepIntExtStart = new INT [IrReps*2];
	if ( mrmos.getIrRepIntExtStartP() )
		memcpy(MOIrRepIntExtStart, mrmos.getIrRepIntExtStartP(),
			2*IrReps*sizeof(INT));

	nInactiveMOs = mrmos.getNumberOfInactiveMOs();
	if ( inactiveMO )
		delete[] inactiveMO;
	inactiveMO = new MOType[nInactiveMOs];
	memcpy(inactiveMO, mrmos.getInactiveMOP(), nInactiveMOs*sizeof(MOType));

	return *this;
}


MRMOs::MRMOs(const MOIrReps &irreps)
{
	maxMO = irreps.getMaxMO();
	IrReps = irreps.getNumberOfIrReps();

	if ( inIrRep )
		delete inIrRep;
	inIrRep = new INT [IrReps];
	memcpy(inIrRep, irreps.getInIrRepP(), IrReps*sizeof(INT));

	if ( MOSymmetry )
		delete MOSymmetry;
	MOSymmetry = new INT [maxMO+1];
	memcpy(MOSymmetry, irreps.getMOSymmetryP(), (maxMO+1)*sizeof(IrRep));

	if ( MOinIrRep )
		delete MOinIrRep;
	MOinIrRep = new MOType[maxMO+1];
	memcpy(MOinIrRep, irreps.getNumberInIrRepP(), (maxMO+1)*sizeof(MOType));

	if ( IrRepStart )
		delete IrRepStart;
	IrRepStart = new MOType[IrReps];
	memcpy(IrRepStart, irreps.getStartMOP(), IrReps*sizeof(MOType));

/*	if ( prodTab )
		delete prodTab;
	prodTab = new IrRep[IrReps*IrReps];
	memcpy(prodTab, irreps.getProdTabP(), IrReps*IrReps*sizeof(IrRep));
*/	
	isIntern = new unsigned char [getMaxMO()+1];
	memset(isIntern, 0, (getMaxMO()+1)*sizeof(unsigned char));
	inIrRepInt = NULL;
	inIrRepExt = NULL;
	MOinIrRepIntExt = NULL;
	MOIrRepIntExtStart = NULL;
	nInactiveMOs = 0;
	inactiveMO = NULL;

}

MRMOs::~MRMOs()
{
	if ( isIntern )
		delete[] isIntern;
	if ( inIrRepInt )
		delete[] inIrRepInt;
	if ( inIrRepExt )
		delete[] inIrRepExt;
	if ( MOinIrRepIntExt )
		delete[] MOinIrRepIntExt;
	if ( MOIrRepIntExtStart )
		delete[] MOIrRepIntExtStart;
	if ( inactiveMO )
		delete[] inactiveMO;
}


void	MRMOs::initIntExt()
{
	inIrRepInt = new INT [IrReps];
	memset(inIrRepInt, 0, IrReps*sizeof(INT));

	inIrRepExt = new INT [IrReps];
	memset(inIrRepExt, 0, IrReps*sizeof(INT));

	for ( MOType i=1 ; i<=maxMO ; i++ )
		if ( isInternal(i) )
			inIrRepInt[getIrRep(i)]++;
		else
			inIrRepExt[getIrRep(i)]++;


	MOinIrRepIntExt = new MOType[maxMO+1];
	MOIrRepIntExtStart = new INT [IrReps*2];
	
INT sum = 0;
	for ( INT i=0 ; i<IrReps ; i++ )
	{	for ( INT j=0 ; j<inIrRep[i] ; j++ )
			MOinIrRepIntExt[j+sum] = 
				j - (j>=getInIrRepInt(i) ? getInIrRepInt(i) : 0);
		MOIrRepIntExtStart[i*2] = sum + 1;
		MOIrRepIntExtStart[i*2 + 1] = sum + getInIrRepInt(i) + 1; 
		sum += inIrRep[i];
	}
}

void	MRMOs::setNumberOfInactiveMOs(INT n)
{
	if ( inactiveMO )
		delete[] inactiveMO;
	inactiveMO = NULL;
	if ( n )
		inactiveMO = new MOType[n];
	nInactiveMOs = n;	
}



void	MRMOs::setInactiveMO(INT i, MOType mo)
{	inactiveMO[i] = mo;	}



ostream& operator<<(ostream & s, const MRMOs & mrmos)
{
	s << *((MOIrReps *) &mrmos) << endl;
	
	s << "internal MOs:" << endl;
	for ( INT i=1 ; i<=mrmos.getMaxMO() ; i++ )
		if ( mrmos.isInternal(i) )
			s << moMapping.getReal(i) << " ";
	s << endl;
	
	s << endl;
	s << "external MOs:" << endl;
	for ( INT i=1 ; i<=mrmos.getMaxMO() ; i++ )
		if ( mrmos.isExternal(i) )
			s << moMapping.getReal(i) << " ";
	s << endl;
	s << endl;
	s << "internal MOs per symmetry:" << endl;
	for ( INT i=0 ; i<mrmos.getNumberOfIrReps() ; i++ )
		s << mrmos.getInIrRepInt(i) << " ";
	s << endl;
INT	nInt = 0;
	for ( INT i=0 ; i<mrmos.getNumberOfIrReps() ; i++ )
		nInt += mrmos.getInIrRepInt(i);
		
	s << "total internal MOs: " << nInt 
		<< " (" << 100.0*nInt/mrmos.maxMO << "%)" << endl;
	s << endl;
	s << "external MOs per symmetry:" << endl;
	for ( INT i=0 ; i<mrmos.getNumberOfIrReps() ; i++ )
		s << mrmos.getInIrRepExt(i) << " ";
	s << endl;

	s << "inactive MOs: ";
	for ( INT i=0 ; i<mrmos.getNumberOfInactiveMOs() ; i++ )
		s << mrmos.getInactiveMO(i) << " ";
	s << endl;
	
	return s;
}
