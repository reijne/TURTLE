//***********************************************************************
//
//	Name:			TupelStructureBase.h
//
//	Description:	
//
//	Author:			Michael Hanrath
//
//	Version:		0.0
//
//	Date:			09.03.1997
//
//
//
//
//
//***********************************************************************


#ifndef __TupelStructureBase_H
#define __TupelStructureBase_H

#include "../../../../config.h"



#include "../../Configuration/Configuration.h"


class TupelStructureBase : 
	public Configuration<MOType>{
public:
	TupelStructureBase() {}
	TupelStructureBase(istream &s);
	TupelStructureBase(IrRep irrep, 
		INT NumberOfOpenShells, INT NumberOfClosedShells, MOType *pShells,
		INT ReferenceFlag = 0);

	TupelStructureBase(IrRep irrep, 
		Configuration<MOType> conf,
		INT ReferenceFlag = 0);
		
	virtual ~TupelStructureBase();
	
	//----------------------------------------------------------------

	IrRep	getIrRep() const;


	INT	isReference() const;

	//----------------------------------------------------------------	

	void	writeToStream(ostream & s) const;

	friend ostream& operator<<(ostream & s, const TupelStructureBase &);
	friend istream& operator>>(istream & s, TupelStructureBase &);

	//----------------------------------------------------------------	

protected:
IrRep	irrep;
INT	ReferenceFlag;		// flag if configuration is reference
};


inline
IrRep	TupelStructureBase::getIrRep() const
{	return irrep;	}


inline
INT	TupelStructureBase::isReference() const
{	return ReferenceFlag;	}



#endif
