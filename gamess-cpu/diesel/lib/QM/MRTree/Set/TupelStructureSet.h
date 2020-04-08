//***********************************************************************
//
//	Name:			TupelStructureSet.h
//
//	Description:	
//
//	Author:			Michael Hanrath
//
//	Version:		0.1
//
//	Date:			30.05.1997
//
//
//
//
//***********************************************************************


#ifndef __TupelStructureSet_H
#define __TupelStructureSet_H

#include "../../../../config.h"



#include "../Base/TupelStructureBase.h"
#include "../../Configuration/Configuration.h"
//#include "extMOsSet.h"
#include "../Container/SetContainerTree.h"


#include "../Base/MRTreeBase.h"


class	extMOsSet;
class	InternalConfsSet;

class TupelStructureSet : 
	public TupelStructureBase, 
	public SetContainerTree<InternalConfsSet, extMOsSet>,
	public MRTreeBase<InternalConfsSet, extMOsSet> {
public:
	TupelStructureSet(InternalConfsSet *parent);
	TupelStructureSet(istream &s);

	TupelStructureSet(IrRep irrep, 
		const Configuration<MOType> & conf,
		InternalConfsSet *parent, INT ReferenceFlag = 0);

	~TupelStructureSet();

	// "virtual" constructor:
//	virtual TupelStructureSet * new_Set();

	void	copy(TupelStructureSet *, INT deep = 0);
	virtual SetContainerTree<InternalConfsSet, extMOsSet> *	clone(INT deep = 0);

	//--------------------------------------------------------------------------
	
	void	add(const Configuration<MOType> &external);
	
	//--------------------------------------------------------------------------

	//	to avoid ambiguity
	TupelStructureSet &	operator &= (TupelStructureSet& b);

	//--------------------------------------------------------------------------

	void	writeToStream(ostream & s) const;

	friend ostream& operator<<(ostream & s, const TupelStructureSet &);
	friend istream& operator>>(istream & s, TupelStructureSet &);

	//----------------------------------------------------------------	

private:
INT	SAFStart;			// start of saf in CI-vector
INT	SAFInc;				// number of SAFs
};

/*
inline 
TupelStructureSet *	TupelStructureSet::
	new_Set()
{	return new TupelStructureSet(NULL);	}
*/

inline 
void	TupelStructureSet::copy(TupelStructureSet *c, INT deep)
{
//	*this = *c;
	(TupelStructureBase &) *this = *c;
	if ( deep )
		SetContainerTree<InternalConfsSet, extMOsSet>::clone(1);
}


inline 
SetContainerTree<InternalConfsSet, extMOsSet> *	TupelStructureSet::clone(INT deep)
{	
	TupelStructureSet * c = new TupelStructureSet(NULL);
	c->copy(this, deep);
	return c;
}


inline
INT operator == (const TupelStructureSet &a, const TupelStructureSet &b)
{	return ((Configuration<MOType>) a) == ((Configuration<MOType>) b);	}

inline
INT operator != (const TupelStructureSet &a, const TupelStructureSet &b)
{	return ((Configuration<MOType>) a) != ((Configuration<MOType>) b);	}

inline
INT operator <= (const TupelStructureSet &a, const TupelStructureSet &b)
{	return ((Configuration<MOType>) a) <= ((Configuration<MOType>) b);	}

inline
TupelStructureSet &	TupelStructureSet::operator &= (TupelStructureSet& b)
{
	*((SetContainerTree<InternalConfsSet, extMOsSet> *) this) &= b;
	return *this;
}


#endif
