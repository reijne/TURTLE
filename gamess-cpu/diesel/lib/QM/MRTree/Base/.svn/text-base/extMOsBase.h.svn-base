//***********************************************************************
//
//	Name:			extMOsBase.h
//
//	Description:	
//
//	Author:			Michael Hanrath
//
//	Version:		0.0
//
//	Date:			11.03.1997
//
//
//
//
//
//***********************************************************************

#ifndef __extMOsBase_H
#define __extMOsBase_H

#include "../../../../config.h"

#include "../../MO/MOType.h"
#include "extEntry.h"
#include "../../Configuration/Configuration.h"
#include "../../Configuration/ConfigurationSAFNr.h"
#include "../Container/Tree.h"
#include "../Container/ContainerIterator.h"
#include "../Container/ContainerTree.h"


#include <iostream>


class	TupelStructureBase;


class extMOsBase :
	public virtual ContainerTree<TupelStructureBase, extEntry> {
public:
	extMOsBase();
	extMOsBase(istream &s);
	extMOsBase(
		INT Syms,
		INT open, 
		INT closed);
		
	virtual ~extMOsBase();

	extMOsBase(const extMOsBase &);
	extMOsBase & operator = (const extMOsBase &);
	//----------------------------------------------------------------

	virtual ConfigurationSAFNr<MOType>	
		getConfigurationSAFNr(const ContainerIterator &) const = 0;


	virtual const extEntry *	getExtEntry(Configuration<MOType> external) const;

	//----------------------------------------------------------------

	
	INT	getNumberOfOpenMOs() const;
	INT	getNumberOfClosedMOs() const;
	INT	getNumberOfTotalMOs() const;
	
	INT	getInSymN(INT sym) const;
	void	setInSymN(INT sym, INT n);

	MOType	getMaxMO() const;
	INT	getNumberOfSymmetries() const;

	INT getSAFStart() const;
	INT getSAFStart(INT number) const;
	INT getSAFInc() const;
	
	//----------------------------------------------------------------

	INT operator == (const extMOsBase &) const;
	INT operator != (const extMOsBase &) const;
	INT operator <= (const extMOsBase &) const;

	//----------------------------------------------------------------	

	void	writeToStream(ostream & s) const;

	friend ostream& operator<<(ostream & s, const extMOsBase &);
	friend istream& operator>>(istream & s, extMOsBase &);

//--------------------------------------------------------------------------

protected:
INT	total;				// sum of open and closed shells
INT	open;				// # of open shells

INT	SAFStart;			// start of SAF in CI-vector
INT	SAFInc;				// number of SAFs

INT	Syms;				// # of symmetries
INT	*inSymN;			// # of configurations in each symmetry

INT	maxMO;				// maximum used MO number
};



inline
INT	extMOsBase::getNumberOfSymmetries() const
{	return	Syms;	}

inline
void	extMOsBase::setInSymN(INT sym, INT n)
{	inSymN[sym] = n;	}

inline
INT	extMOsBase::getInSymN(INT sym) const
{	return inSymN[sym];	}

inline
INT	extMOsBase::getNumberOfOpenMOs() const
{	return open;	}

inline
INT	extMOsBase::getNumberOfClosedMOs() const
{	return total-open;	}

inline
INT	extMOsBase::getNumberOfTotalMOs() const
{	return total;	}

inline
INT extMOsBase::getSAFStart() const
{	return SAFStart;	}

inline
INT extMOsBase::getSAFStart(INT number) const
{	return SAFStart + number*SAFInc;	}

inline
INT extMOsBase::getSAFInc() const
{	return SAFInc;	}

inline
MOType	extMOsBase::getMaxMO() const
{	return maxMO;	}



inline
INT extMOsBase::operator == (const extMOsBase &a) const
{
	return total==a.getNumberOfTotalMOs() &&
		open==a.getNumberOfOpenMOs();
}

inline
INT extMOsBase::operator != (const extMOsBase &a) const
{
	return total!=a.getNumberOfTotalMOs() ||
		open!=a.getNumberOfOpenMOs();
}

inline
INT extMOsBase::operator <= (const extMOsBase &a) const
{
	if ( total<a.getNumberOfTotalMOs() )
		return 1;
	if ( open<=a.getNumberOfOpenMOs() )
		return 1;
	return 0;
}



#endif
