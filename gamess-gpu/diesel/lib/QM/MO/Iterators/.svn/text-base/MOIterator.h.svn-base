//***********************************************************************
//
//	Name:			MOIterator.h
//
//	Description:	MOs classified by internal/external status and
//					irrep to be used as creators or annihilators on
//					internal configuration rests
//
//	Author:			Michael Hanrath
//
//	Version:		0.0
//
//	Date:			20.09.1998
//
//
//
//
//
//***********************************************************************


#ifndef __MOIterator_h
#define __MOIterator_h

#include "../../../../config.h"


#include "../../Symmetry/IrRep.h"
#include "../MOType.h"
#include "../MRMOs.h"
#include "../../Configuration/Configuration.h"


/* FD class ostream; */



class MOIterator {
public:
	//	all n-tupels
	MOIterator(
		INT nTupel, 
		const MRMOs *mrmos, 
		INT isExternal);

	//	all n-tupels with specific product irrep
	MOIterator(
		INT nTupel, 
		const MRMOs *mrmos, 
		INT isExternal, 
		IrRep productIrRep);


	~MOIterator();
	
	operator Configuration<MOType>() const;
	
	
	void	next();
	INT isEnd() const;
	
	INT	changedCase() const;
	
	MOType	getMO(INT) const;
	MOType	getMO(MOType, INT) const;
	MOType	getMOEnd(INT) const;
	
	void	skip(INT);
	
	INT	getNumberOfOpenShells() const;
	
	friend ostream & operator << (ostream &, const MOIterator &);
	
	
private:
	MOIterator(const MOIterator &);
	MOIterator & operator = (const MOIterator &);

	void	init();
	
	INT	irrepIter();
	INT	moIter();
	
	
INT	nTupel;						// n-tupel
const MRMOs	*mrmos;				// MRMOs
INT	isExternal;					// flag if MOs are from isExternal set
IrRep	productIrRep;			// product irrep
INT	specificIrrep;				// flag if specific irrep is ordered


INT	end;						// flag if no further iteration
INT	changed;					// flag if interaction case has changed 
								// (for external cases only)
								//   - tupels from different irrep
								//   - different tupel structure


//===========================================================================
// iterator states:

//---------------------------------------------------------------------------
// running MOs:
MOType	*mos;
INT	iMOs;

//---------------------------------------------------------------------------
// irreps:
INT	lIrRep;						// internal loop counter for irreps
IrRep	*irrep;					// internal status of irreps ("upper triangle")
};


inline
INT	MOIterator::isEnd() const
{	return end;	}

inline
INT	MOIterator::changedCase() const
{	return changed;	}

inline
INT	MOIterator::getNumberOfOpenShells() const
{
INT	n = 1;
	for ( INT i=1 ; i<nTupel ; i++ )
		if ( mos[i]==mos[i-1] )
			n--;
		else
			n++;
	return n;

}

inline
MOType	MOIterator::getMO(INT i) const
{	return	mos[i];	}

inline
MOType	MOIterator::getMO(MOType mo, INT iMOs) const
{
	if ( iMOs>0 && irrep[iMOs-1]==irrep[iMOs] )
		return mo + (iMOs>2 ? irrep[iMOs-2]==irrep[iMOs] : 0);
	else
		return mrmos->getIrRepIntExtStart(irrep[iMOs], isExternal);
}


inline
MOType	MOIterator::getMOEnd(INT i) const
{	return	mrmos->getIrRepIntExtEnd(irrep[i], isExternal);	}


#endif
