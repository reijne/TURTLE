//***********************************************************************
//
//	Name:			InternalConfsDiag.h
//
//	Description:	
//
//	Author:			Michael Hanrath
//
//	Version:		0.0
//
//	Date:			21.06.1996
//
//
//
//
//
//***********************************************************************


#ifndef __InternalConfsDiag_H
#define __InternalConfsDiag_H

#include "../../../../config.h"


#include "../Base/InternalConfsBase.h"

#include "../Diag/TupelStructureDiag.h"
#include "../Diag/DiagTreeBase.h"

template <class TMOType> class Configuration;
class NExternalsDiag;
class InternalConfsSet;


class InternalConfsDiag : 
	public InternalConfsBase,
	public DiagTreeBase<NExternalsDiag, TupelStructureDiag> {
public:
	InternalConfsDiag(INT nExt, INT n, NExternalsDiag *parent) :
                Tree<NExternalsDiag>(parent),
		IndexedContainer<TupelStructureDiag>(n),
                InternalConfsBase(nExt) {}
	InternalConfsDiag(istream &s);

	InternalConfsDiag(const InternalConfsDiag &);

	InternalConfsDiag(const InternalConfsSet &, NExternalsDiag *, INT referenceFlag);
	
	~InternalConfsDiag() {}
	
	//----------------------------------------------------------------	

	void	writeToStream(ostream & s) const;

	friend ostream& operator<<(ostream & s, const InternalConfsDiag &);
	friend istream& operator>>(istream & s, InternalConfsDiag &);

	//----------------------------------------------------------------	


private:
};


#endif
