//***********************************************************************
//
//	Name:			MOList.h
//
//	Description:	MOs classified by internal/external status and
//					irrep to be used as creators or annihilators on
//					internal configuration rests
//
//	Author:			Michael Hanrath
//
//	Version:		0.0
//
//	Date:			05.02.1997
//
//
//
//
//
//***********************************************************************

#ifndef __MOList_h
#define __MOList_h

#include "../../../../config.h"

#include <string>
#include <iostream>

#include "../../Symmetry/IrRep.h"
#include "../MOType.h"
#include "../MRMOs.h"

class MOList {
public:
	MOList(INT nTupel, const MRMOs *mrmos, INT *isInternal, IrRep *irrep);
	~MOList();


	INT	getNTupel() const;
	INT getIsInternal(INT i) const;
	IrRep getIrRep(INT i) const;
	IrRep getIrRepExt(INT i) const;


	INT	getNumberOfExternalIrReps() const;
	IrRep	getExternalIrRep(INT i) const;
	INT	getInExternalIrRep(INT i) const;


	MOType getIrRepExtStart(INT i) const;
	MOType getIrRepExtEnd(INT i) const;
	
	const MRMOs	*getMRMOs() const;

	INT	getNumberOfInternals() const;
	INT	getNumberOfExternals() const;
		
	INT	getTotalNumber() const;

	//----------------------------------------------------------------	

	friend ostream& operator<<(ostream & s, const MOList &);

	//----------------------------------------------------------------

private:
INT	nTupel;				//	number of MOs in a tupel (e.g. r=1, rs=2, rst=3, ...)
INT	nInternal;			//	number of internal MOs
IrRep	*irrep;			//	irrep of MOs, internal irreps first
const MRMOs	*mrmos;		//	pointer to information about intern/external
						//	structure and irreps
						
INT	nExtIrReps;			//	number of external irreps
IrRep	*extIrRep;		//	nth external irrep
INT	*inExtIrRep;		//	number of external electrons in specific irrep

};


inline
MOList::MOList(INT _nTupel, const MRMOs *_mrmos, INT *_isInternal, IrRep *_irrep)
{
	nTupel = _nTupel;
	mrmos = _mrmos;
	nInternal = 0;
	irrep = new IrRep[nTupel];
	for ( INT i=0 ; i<nTupel ; i++ )
		if ( _isInternal[i] )
			irrep[nInternal++] = _irrep[i];


	nExtIrReps = 0;
	extIrRep = new IrRep[mrmos->getNumberOfIrReps()];	
	inExtIrRep = new INT[mrmos->getNumberOfIrReps()];	

	
IrRep	prev = -1;
INT	j = nInternal;
	for ( INT i=0 ; i<nTupel ; i++ )
		if ( !_isInternal[i] )
		{
			irrep[j++] = _irrep[i];
			if ( prev != _irrep[i] )
			{
				inExtIrRep[nExtIrReps] = 1;
				extIrRep[nExtIrReps++] = prev = _irrep[i];
			}
			else
				inExtIrRep[nExtIrReps-1]++;
		}
}

inline
MOList::~MOList()
{
	if ( nTupel )
		delete[] irrep;
		
	delete[] extIrRep;
	delete[] inExtIrRep;
}

inline
INT	MOList::getNTupel() const
{	return nTupel;	}

inline
INT MOList::getIsInternal(INT i) const
{	return (i<nInternal);	}

inline
IrRep MOList::getIrRep(INT i) const
{	return irrep[i];	}	

inline
IrRep MOList::getIrRepExt(INT i) const
{	return irrep[nInternal + i];	}	

inline
const MRMOs	*MOList::getMRMOs() const
{	return mrmos;	}

inline
MOType MOList::getIrRepExtStart(INT i) const
{	return mrmos->getIrRepIntExtStart(irrep[nInternal + i], 1);	}

inline
MOType MOList::getIrRepExtEnd(INT i) const
{	return mrmos->getIrRepIntExtEnd(irrep[nInternal + i], 1);	}

inline
INT	MOList::getNumberOfInternals() const
{	return nInternal;	}

inline
INT	MOList::getNumberOfExternals() const
{	return nTupel - nInternal;	}

inline
INT	MOList::getNumberOfExternalIrReps() const
{	return	nExtIrReps;	}

inline
IrRep	MOList::getExternalIrRep(INT i) const
{	return	extIrRep[i];	}

inline
INT	MOList::getInExternalIrRep(INT i) const
{	return	inExtIrRep[i];	}



#endif
