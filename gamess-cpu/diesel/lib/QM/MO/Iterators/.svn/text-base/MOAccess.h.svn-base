//***********************************************************************
//
//	Name:			MOAccess.h
//
//	Description:	address lists for fast access
//
//	Author:			Michael Hanrath
//
//	Version:		0.0
//
//	Date:			26.09.1998
//
//
//
//
//
//***********************************************************************


#ifndef __MOAccess_h
#define __MOAccess_h

#include "../../../../config.h"


#include "../MOType.h"


class MRMOs;

class NExternalsDiag;

class MOAccess {
public:
	MOAccess(INT nTupel, const MRMOs *mrmos, INT maxInternalOpen,
	const NExternalsDiag &);
	~MOAccess();
	
	INT	getConfOffset(MOType mo) const;
	INT	getConfOffset(MOType mo1, MOType mo2) const;

	INT	getSAFOffset(INT internalOpen, MOType mo) const;
	INT	getSAFOffset(INT internalOpen, MOType mo1, MOType mo2) const;
//	INT	operator () (INT open, const MOType *mo) const;

	
private:
INT nTupel;
INT	maxMO;
INT maxInternalOpen;
INT	*pConf;
INT	**pSAF;
const NExternalsDiag & Psi0;				// configuration tree to be excluded
};



inline
INT	MOAccess::getConfOffset(MOType mo) const
{	return pConf[mo-1];	}


inline
INT	MOAccess::getConfOffset(MOType mo1, MOType mo2) const
{
	if ( mo2>mo1 )
		return pConf[(mo2-1)*maxMO + mo1 - 1];
	else
		return pConf[(mo1-1)*maxMO + mo2 - 1];
}


inline
INT	MOAccess::getSAFOffset(INT internalOpen, MOType mo) const
{	return pSAF[internalOpen][mo-1];	}


inline
INT	MOAccess::getSAFOffset(INT internalOpen, MOType mo1, MOType mo2) const
{
	if ( mo2>mo1 )
		return pSAF[internalOpen][(mo2-1)*maxMO + mo1 - 1];
	else
		return pSAF[internalOpen][(mo1-1)*maxMO + mo2 - 1];
}

/*
inline
INT	MOAccess::operator () (INT open, const MOType *mo) const
{

}
*/


#endif
