//***********************************************************************
//
//	Name:			TupelStructureSel.h
//
//	Description:	
//
//	Author:			Michael Hanrath
//
//	Version:		0.1
//
//	Date:			21.06.1996
//					09.03.1997
//
//
//
//
//
//***********************************************************************


#ifndef __TupelStructureSel_H
#define __TupelStructureSel_H

#include "../../../../config.h"



#include "../Base/TupelStructureBase.h"
#include "../../Configuration/Configuration.h"
#include "extMOsSel.h"
#include "../Container/IndexedContainerTree.h"


#include "../Base/MRTreeBase.h"


class	extMOsSel;
class	InternalConfsSel;
class	MOEquivalence;

class TupelStructureSel : 
	public TupelStructureBase, 
	public IndexedContainerTree<InternalConfsSel, extMOsSel>,
	public MRTreeBase<InternalConfsSel, extMOsSel> {
public:
	TupelStructureSel(InternalConfsSel *parent);
	TupelStructureSel(istream &s);

	TupelStructureSel(INT nElectrons, IrRep irrep, 
		INT NumberOfOpenShells, INT NumberOfClosedShells, MOType *pShells,
		InternalConfsSel *parent, INT ReferenceFlag = 0);

	TupelStructureSel(INT number, IrRep irrep, 
		Configuration<MOType> conf,
		InternalConfsSel *parent, INT ReferenceFlag = 0);

	~TupelStructureSel();

	// "virtual" constructor:
//	virtual TupelStructureSel * new_Set();
	

	void	copy(TupelStructureSel *, INT deep = 0);
	virtual TupelStructureSel *	clone(INT deep = 0);
	
	
	//----------------------------------------------------------------	

	void	addExtEntry(const Configuration<MOType> & external,
		const RootEnergies &);
	const extEntry *	getExtEntry(Configuration<MOType> external) const;

	//--------------------------------------------------------------------------

	void	threshCut(EnergyType Threshold, INT divideByCSFs);

	//--------------------------------------------------------------------------

	friend INT operator == (const TupelStructureSel &,
		const TupelStructureSel &);
	friend INT operator != (const TupelStructureSel &,
		const TupelStructureSel &);
	friend INT operator <= (const TupelStructureSel &,
		const TupelStructureSel &);


	INT	operator |= (TupelStructureSel) { return 0; }

	//----------------------------------------------------------------	

/*
	void	add(
		const Configuration<MOType> & conf,
		INT	roots,
		const EnergyType *energy);
*/

	//----------------------------------------------------------------	

	INT getSAFStart() const;
	INT getSAFInc() const;
	
	INT	init(INT, const IrRep *);
		
	//----------------------------------------------------------------	

	void	writeToStream(ostream & s) const;

	friend ostream& operator<<(ostream & s, const TupelStructureSel &);
	friend istream& operator>>(istream & s, TupelStructureSel &);

	//----------------------------------------------------------------	

private:
INT	SAFStart;			// start of saf in CI-vector
INT	SAFInc;				// number of SAFs
};


/*
inline 
TupelStructureSel *	TupelStructureSel::
	new_Set()
{	return new TupelStructureSel(NULL);	}
*/

inline 
void	TupelStructureSel::copy(TupelStructureSel *c, INT deep)
{
	if ( deep )
	{
		*this = *c;
		// !!!!!!!!!! deep copy missing
//		IndexedContainerTree<InternalConfsSel, extMOsSel>::clone(1);
	}
	else
		*this = *c;

}


inline 
TupelStructureSel *	TupelStructureSel::clone(INT deep)
{	
	TupelStructureSel * c = new TupelStructureSel(NULL);
	c->copy(this, deep);
	return c;
}


inline
INT TupelStructureSel::getSAFStart() const
{	return SAFStart;	}

inline
INT TupelStructureSel::getSAFInc() const
{	return SAFInc;	}

inline
INT operator == (const TupelStructureSel &a, const TupelStructureSel &b)
{	return ((Configuration<MOType>) a) == ((Configuration<MOType>) b);	}

inline
INT operator != (const TupelStructureSel &a, const TupelStructureSel &b)
{	return ((Configuration<MOType>) a) != ((Configuration<MOType>) b);	}

inline
INT operator <= (const TupelStructureSel &a, const TupelStructureSel &b)
{	return ((Configuration<MOType>) a) <= ((Configuration<MOType>) b);	}


#endif
