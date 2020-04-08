//***********************************************************************
//
//	Name:			MRMOs.h
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

#ifndef __MRMOS_H
#define __MRMOS_H

#include "../../../config.h"

#include "MOType.h"
#include "MOIrReps.h"

#include <iostream>
using std::ostream;

class MRMOs : public MOIrReps {
public:
	MRMOs();
	MRMOs(Fort31File f31);
	~MRMOs();
	
	// copy-constructor
	MRMOs(const MRMOs &);

	// Zuweisungs-Operator
	MRMOs & operator = (const MRMOs &);

	// type conversion
	MRMOs(const MOIrReps &);

	//----------------------------------------------------------------	

	IrRep	operator [] (INT);

	//----------------------------------------------------------------	

	INT	isInternal(MOType mo) const;
	INT	isExternal(MOType mo) const;
	
	void	setInternal(MOType mo, INT flag = 1);
	void	setExternal(MOType mo, INT flag = 1);

	INT	getInIrRepInt(IrRep i) const;
	INT	getInIrRepExt(IrRep i) const;

	MOType	getNumberInIrRepIntExt(MOType mo) const;

	MOType	getIrRepIntExtStart(IrRep, INT ext) const;
	MOType	getIrRepIntExtEnd(IrRep, INT ext) const;
	INT	getNumberOfIrRepIntExt(IrRep, INT ext) const;
	
	INT	getNumberOfInactiveMOs() const;
	MOType	getInactiveMO(INT i) const;

	void	setNumberOfInactiveMOs(INT n);
	void	setInactiveMO(INT i, MOType mo);
	
	const MOType *getInactiveMOP() const;
	
	//----------------------------------------------------------------	

	const unsigned char *	getInternP() const;
	const INT *	getInIrRepIntP() const;
	const INT *	getInIrRepExtP() const;
	const MOType *	getNumberInIrRepIntExtP() const;
	const INT *	getIrRepIntExtStartP() const;

	//----------------------------------------------------------------	

	void	initIntExt();

	//----------------------------------------------------------------	

	friend ostream& operator<<(ostream & s, const MRMOs &);

	//----------------------------------------------------------------
		
protected:
unsigned char *isIntern;	// flag if MO is internal
INT	*inIrRepInt;			// number of internal MOs in irrep
INT	*inIrRepExt;			// number of external MOs in irrep
MOType	*MOinIrRepIntExt;	// translation table MO -> INT/ext MO in irrep
INT	*MOIrRepIntExtStart;	// translation table MO -> INT/ext start in irrep
INT	nInactiveMOs;			// number of inactive MOs
MOType	*inactiveMO;			// list of inactive MOs
};



inline
INT	MRMOs::isInternal(MOType mo) const
{	return	isIntern[mo];	}

inline
INT	MRMOs::isExternal(MOType mo) const
{	return	!isIntern[mo];	}

inline
void	MRMOs::setInternal(MOType mo, INT flag)
{	isIntern[mo] = flag;	}

inline
void	MRMOs::setExternal(MOType mo, INT flag)
{	isIntern[mo] = !flag;	}

inline
INT	MRMOs::getInIrRepInt(IrRep i) const
{	return inIrRepInt[i];	}

inline
INT	MRMOs::getInIrRepExt(IrRep i) const
{	return inIrRepExt[i];	}

inline
const unsigned char *	MRMOs::getInternP() const
{	return isIntern;	}

inline
const INT *	MRMOs::getInIrRepIntP() const
{	return inIrRepInt;	}

inline
const INT *	MRMOs::getInIrRepExtP() const
{	return inIrRepExt;	}

inline
const MOType *	MRMOs::getNumberInIrRepIntExtP() const
{	return MOinIrRepIntExt; }

inline
const INT *	MRMOs::getIrRepIntExtStartP() const
{	return MOIrRepIntExtStart;}

inline
MOType	MRMOs::getNumberInIrRepIntExt(MOType mo) const
{	return MOinIrRepIntExt[mo];	}


inline
MOType	MRMOs::getIrRepIntExtStart(IrRep i, INT ext) const
{
	if ( ext )
		return MOIrRepIntExtStart[i*2+1];
	else
		return MOIrRepIntExtStart[i*2+0];
}


inline
MOType	MRMOs::getIrRepIntExtEnd(IrRep i, INT ext) const
{
	if ( ext )
	{	if ( i+1>=IrReps )
			return maxMO;
		else
			return MOIrRepIntExtStart[i*2+2] - 1;
	}
	else
		return MOIrRepIntExtStart[i*2+1] - 1;
}

inline
INT	MRMOs::getNumberOfIrRepIntExt(IrRep i, INT ext) const
{
	if ( ext )
	{	if ( i+1>=IrReps )
			return maxMO - MOIrRepIntExtStart[i*2+1] + 1;
		else
			return MOIrRepIntExtStart[i*2+2] - MOIrRepIntExtStart[i*2+1];
	}
	else
		return MOIrRepIntExtStart[i*2+1] - MOIrRepIntExtStart[i*2];
}



inline
INT	MRMOs::getNumberOfInactiveMOs() const
{	return nInactiveMOs;	}



inline
MOType	MRMOs::getInactiveMO(INT i) const
{	return	inactiveMO[i];	}


inline
const MOType *MRMOs::getInactiveMOP() const
{	return	inactiveMO;	}




#endif
