//***********************************************************************
//
//	Name:			MatchingExternalMOIterator.cc
//
//	Description:	
//
//	Author:			Michael Hanrath
//
//	Version:		0.0
//
//	Date:			24.08.1998
//
//
//
//
//
//***********************************************************************




#include "MatchingExternalMOIterator.h"


MatchingExternalMOIterator::MatchingExternalMOIterator(
	const extMOsDiag *_MOA,
	const extMOsDiag *_MOB,
	INT _maxExcitationPossible,
	INT _maxExcitationAllowed,
	INT isDiag,
	INT	_la
	)
{
	MOA = _MOA;
	MOB = _MOB;
	la = _la;
	maxExcitationPossible = _maxExcitationPossible;
	maxExcitationAllowed = _maxExcitationAllowed;
	if ( isDiag )
	{
		ind = la-1;
	}
	else
	{
		ind = -1;
	}
	maxInd = MOB->getNumberOfElements();

/*	cout << "ind=" << ind << ", maxInd=" << maxInd << endl;
	cout << "maxEcitationPossible = " << maxExcitationPossible << endl;
	cout << "maxEcitationAllowed = " << maxExcitationAllowed << endl;
*/
	// no match possible
/*	if ( maxExcitationPossible - 
		(2*MOB->getNumberOfClosedMOs() + MOB->getNumberOfOpenMOs()) >
		maxExcitationAllowed )
	{
		maxInd = 0;
		maxExcitationPossible = -1;
	}
	else
*/
	if ( maxExcitationPossible>maxExcitationAllowed )
	{


		maskEnd = 
			((MOA->getNumberOfTotalMOs()==2) << 1) +
			((MOB->getNumberOfTotalMOs()==2) << 0) + 1;
			
		maskInc = 1 + 
			(MOA->getNumberOfTotalMOs()==2 && MOB->getNumberOfTotalMOs()==1);

		mask = -maskInc;
/*		cout << "maskEnd = " << maskEnd << endl;
		cout << "maskInc = " << maskInc << endl;
		cout << "MOA=" << MOA->getNumberOfTotalMOs() << endl << *MOA << endl;
		cout << "la=" << la << endl;
		cout << "MOB=" << MOB->getNumberOfTotalMOs() << endl << *MOB << endl;
*/

		lb = 0;
	}
}




