//***********************************************************************
//
//	Name:			StoneyToPlainIndex.h
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



#ifndef __STONEYTOPLAININDEX_H
#define __STONEYTOPLAININDEX_H

#include "config.h"


#include "PlainIndex.h"

#include "../IntegralIndex/TwoElectronIntegralIndex.h"
#include "../MO/MOIrReps.h"




class StoneyToPlainIndex : public PlainIndex, private MOIrReps {
public:
	StoneyToPlainIndex(const char *Fort31FileName);
	~StoneyToPlainIndex();


	//----------------------------------------------------------------

	INT	getNumberOfIntegrals() const;

	void	set(const TwoElectronIntegralIndex<MOType> &);

	//----------------------------------------------------------------

	friend ostream& operator<<(ostream & s, const StoneyToPlainIndex &);

	//----------------------------------------------------------------
private:
	MOType	getNumberInIrRep(MOType mo) const;
	MOType	getNumberInIrRep(MOType mo, IrRep i) const;

const INT	maxIrRepBlocks = 666;
INT	*irRepBlockStart;
INT	*inIrRepSum;
INT	integrals;
};


inline
MOType	StoneyToPlainIndex::getNumberInIrRep(MOType mo) const
{	return	mo - inIrRepSum[getIrRep(mo)];	}

inline
MOType	StoneyToPlainIndex::getNumberInIrRep(MOType mo, IrRep i) const
{	return	mo - inIrRepSum[i];	}





#endif
