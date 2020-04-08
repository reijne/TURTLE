//***********************************************************************
//
//	Name:			MatchingExternalMOIterator.h
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


#ifndef __MatchingExternalMOIterator_h
#define __MatchingExternalMOIterator_h

#include "../../../../config.h"


#include "../../MO/MOType.h"


class extMOsDiag;

class MatchingExternalMOIterator {
public:
	MatchingExternalMOIterator(
		const extMOsDiag *MOA,
		const extMOsDiag *MOB,
		INT maxExcitationPossible,
		INT maxExcitationAllowed,
		INT isDiag,
		INT	la
		);
	~MatchingExternalMOIterator() {};

	INT	next();
	INT	getIndex() const;
	
	
private:
const extMOsDiag *MOA;
const extMOsDiag *MOB;
INT	la;
INT	lb;
INT	maxExcitationPossible;
INT	maxExcitationAllowed;
INT	ind;
INT	maxInd;
INT	mask;
INT	maskEnd;
INT	maskInc;
MOType	mo1, mo2;
};


#include "extMOsDiag.h"


inline
INT	MatchingExternalMOIterator::next()
{
	if ( maxExcitationPossible<=maxExcitationAllowed )
		return ++ind<maxInd;
	else
	{
		if ( mask>=0 )
		{
//			cout << "lb (++) = " << lb << endl;

			while ( ++lb<maxInd && MOB->getOpenMO(lb, mask & 1, mask & 1) == mo1 )
			{
//				cout << "mask=" << mask << ", lb=" << lb << ", " << mo1 << endl;
//				cout << MOA->getOpenMO(la, !(mask & 2), 0) << ":::" <<
//					MOB->getOpenMO(lb, !(mask & 1), !!(mask & 1)) << endl;

				// avoid double generation in rs|tu case with r=t and s=u
				if ( lb>=0 && mask>=2 && maskInc==1 &&
					MOA->getOpenMO(la, !(mask & 2), 0) == MOB->getOpenMO(lb, !(mask & 1), !!(mask & 1)) )
					;
				else
				{
					ind = MOB->getIndex(lb, mask & 1);
					return 1;
				}
			}
		}
		do
		{
			mask += maskInc;
			if ( mask>=maskEnd )
				return 0;
			lb = MOB->findOpenMO(mask & 1, 
				mo1=MOA->getOpenMO(la, !!(mask & 2)));
				

			// avoid double generation in rs|tu case with r=t and s=u
			if ( lb>=0 && mask>=2 && maskInc==1 &&
				MOA->getOpenMO(la, !(mask & 2), 0) == MOB->getOpenMO(lb, !(mask & 1), !!(mask & 1)) )
			{
				if ( ++lb>=maxInd || MOB->getOpenMO(lb, mask & 1, mask & 1) != mo1 )
					lb = -1;
			}
//			cout << "mask=" << mask << ", lb=" << lb << endl;
//			if ( lb>=0 )
//				cout << MOB->getOpenMO(lb, mask & 1, mask & 1) << " " << mo1 << endl;
		} while ( lb<0 );
		
		ind = MOB->getIndex(lb, mask & 1);
		return 1;
	
	}
}

inline
INT	MatchingExternalMOIterator::getIndex() const
{	return ind;	}




#endif
