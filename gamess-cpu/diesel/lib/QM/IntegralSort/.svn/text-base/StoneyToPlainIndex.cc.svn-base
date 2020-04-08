//***********************************************************************
//
//	Name:			StoneyToPlainIndex.cc
//
//	Description:	index conversion from stoney (ij|kl) to plain index
//
//	Author:			Michael Hanrath
//
//	Version:		0.0
//
//	Date:			30.08.1996
//
//
//
//
//
//***********************************************************************




#include "StoneyToPlainIndex.h"


#include "../IO/Fortran/Fort31SecondRecord.h"
#include "../IntegralIndex/TwoElectronIntegralTriadeIndex.h"



StoneyToPlainIndex::StoneyToPlainIndex(const char *Fort31FileName) :
	MOIrReps(Fort31FileName)
{
Fort31SecondRecord	f31(Fort31FileName);

	irRepBlockStart = new INT[maxIrRepBlocks];
	
	memcpy(irRepBlockStart, &f31.getNIT(0), maxIrRepBlocks*sizeof(INT));

	integrals = f31.getNN();
//	for ( INT i=0 ; i<maxIrRepBlocks ; i++ )
//		cout << irRepBlockStart[i] << endl;

	inIrRepSum = new INT[IrReps];
	
	inIrRepSum[0] = 0;
	for ( INT i=1 ; i<IrReps ; i++ )
		inIrRepSum[i] = inIrRepSum[i-1] + inIrRep[i-1];
}


StoneyToPlainIndex::~StoneyToPlainIndex()
{
	delete inIrRepSum;
	delete irRepBlockStart;
}



INT	StoneyToPlainIndex::getNumberOfIntegrals() const
{	return integrals;	}





void	StoneyToPlainIndex::set(
	const TwoElectronIntegralIndex<MOType> & twoind)
{
	// calculate start address in block
	index = irRepBlockStart[
		iover2(
			iover2(getIrRepP1(twoind.getI())) + 
			getIrRepP1(twoind.getJ())
			) +
		iover2(getIrRepP1(twoind.getK())) + 
		getIrRepP1(twoind.getL())];
		

TwoElectronIntegralTriadeIndex	tt(twoind);
	cout << "vorher:" << tt << endl;
	tt.setTriade();
	cout << "nacher:" << tt << endl;

	// calculate address within block
	if ( getIrRep(twoind.getI())==getIrRep(twoind.getJ()) )
	{	if ( getIrRep(twoind.getI())==getIrRep(twoind.getK()) )
		{	// Fall: (AAAA), all same irrep
			index =	iover2(
							iover2(getNumberInIrRep(tt.getI())) + 
							getNumberInIrRep(tt.getJ())
							) +
						iover2(getNumberInIrRep(tt.getK())) + 
						getNumberInIrRep(tt.getL());
		}
		else
		{	// Fall: (AABB), two different irreps
		}
	}
	else
	{	// Fall: (ABCD)-Irrep, all different irreps
	}
//	index += 
		
	cout << "startindex = " << index << endl;
}



ostream& operator<<(ostream & s, const StoneyToPlainIndex & stind)
{
	s << "number of two electron integrals is " 
		<< stind.getNumberOfIntegrals() << endl;
	
}
