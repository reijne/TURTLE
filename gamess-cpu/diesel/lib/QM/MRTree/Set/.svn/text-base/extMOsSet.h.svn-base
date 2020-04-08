//***********************************************************************
//
//	Name:			extMOsSet.h
//
//	Description:	
//
//	Author:			Michael Hanrath
//
//	Version:		0.0
//
//	Date:			30.05.1997
//
//
//
//
//
//***********************************************************************

#ifndef __extMOsSet_H
#define __extMOsSet_H

#include "../../../../config.h"

#include "../Base/extMOsBase.h"

#include "../../MO/MOType.h"
#include "../../Configuration/Configuration.h"
#include "../../Configuration/ConfigurationSAFNr.h"

//#include "TupelStructureSet.h"
#include "../MRTreeIterator.h"

#include "../Base/extEntry.h"
#include "../Container/SetContainerTree.h"

#include <iostream>

#include "../EnergyType.h"


class	TupelStructureSet;


class extMOsSet  : 
	public SetContainerTree<TupelStructureBase, extEntry>,
	public extMOsBase {
public:
	extMOsSet(TupelStructureBase *parent = NULL):
		Tree<TupelStructureBase>(parent) {}

	extMOsSet(INT open, INT closed, TupelStructureBase *parent = NULL):
		Tree<TupelStructureBase>(parent),
		extMOsBase(1, open, closed) {}

	extMOsSet(istream &s):
                Tree<TupelStructureBase>(NULL),
                SetContainer<extEntry>(s),
                extMOsBase(s) {}
//		extMOsBase(s),
//		SetContainer<extEntry>(s),
//		Tree<TupelStructureBase>(NULL) {}
	
	~extMOsSet();

	// "virtual" constructor:
//	virtual extMOsSet * new_Set();

	void	copy(extMOsSet *, INT deep = 0);
	virtual SetContainerTree<TupelStructureBase, extEntry> *	clone(INT deep = 0);
	
	
	
	//----------------------------------------------------------------
	
	void	add(const Configuration<MOType> &external);
	
	//----------------------------------------------------------------

	ConfigurationSAFNr<MOType>	getConfigurationSAFNr(INT n) const;

	ConfigurationSAFNr<MOType>	// to kill virtual function
		getConfigurationSAFNr(const ContainerIterator &) const;

	ConfigurationSAFNr<MOType>	
		getConfigurationSAFNr(Pix pix) const;

	TupelStructureSet * getParent() const;

	INT calcNumberOfLeaves() const;
	INT getNumberOfLeaves() const;
	
	//--------------------------------------------------------------------------

	INT	getSelectionFlags(const ContainerIterator &n = (INT) 0) const;

	const RootEnergies &	getEnergy(
		const ContainerIterator &n = (INT) 0) const;

	void	setOpenMO(const ContainerIterator &n, INT i, MOType mo);
	void	setClosedMO(const ContainerIterator &n, INT i, MOType mo);
	
	MOType	getOpenMO(const ContainerIterator &n, INT i) const;
	MOType	getClosedMO(const ContainerIterator &n, INT i) const;

	//--------------------------------------------------------------------------

	INT	init(INT, const IrRep *);

	//--------------------------------------------------------------------------

	INT operator == (const extMOsSet &) const;
	INT operator != (const extMOsSet &) const;
	INT operator <= (const extMOsSet &) const;

	//----------------------------------------------------------------
	
	void	writeToStream(ostream & s) const;

	friend ostream& operator<<(ostream & s, const extMOsSet &);
	friend istream& operator>>(istream & s, extMOsSet &);

//--------------------------------------------------------------------------

private:
};


/*
inline 
extMOsSet *	extMOsSet::
	new_Set()
{	return new extMOsSet();	}
*/


inline 
void	extMOsSet::copy(extMOsSet *c, INT deep)
{
//	*this = *c;
	(extMOsBase &) *this = *c;
	if ( deep )
		SetContainerTree<TupelStructureBase, extEntry>::clone(1);
}


inline 
SetContainerTree<TupelStructureBase, extEntry> *	extMOsSet::clone(INT deep)
{	
	extMOsSet * c = new extMOsSet();
	c->copy(this, deep);
	return c;
}

inline
TupelStructureSet * extMOsSet::getParent() const
{	return ((TupelStructureSet *) extMOsBase::getParent());	}


inline
INT	extMOsSet::getNumberOfLeaves() const
{	return length();  }

inline
INT	extMOsSet::calcNumberOfLeaves() const
{	return length();  }

//--------------------------------------------------------------------------

inline
ConfigurationSAFNr<MOType>	extMOsSet::
	getConfigurationSAFNr(const ContainerIterator & iter) const
{	return getConfigurationSAFNr(iter.pix);	}

//--------------------------------------------------------------------------


inline
INT extMOsSet::operator == (const extMOsSet &a) const
{	return extMOsBase::operator == (a);	}

inline
INT extMOsSet::operator != (const extMOsSet &a) const
{	return extMOsBase::operator != (a);	}

inline
INT extMOsSet::operator <= (const extMOsSet &a) const
{	return extMOsBase::operator <= (a);	}

//----------------------------------------------------------------


inline
const RootEnergies &	extMOsSet::getEnergy(const ContainerIterator &n) const
{	return ((extMOsSet *) this)->operator [] (n)->getEnergy();	}

inline
INT	extMOsSet::getSelectionFlags(const ContainerIterator &n ) const
{	return ((extMOsSet *) this)->operator [] (n)->getSelectionFlags();	}

//----------------------------------------------------------------

inline
void	extMOsSet::setOpenMO(const ContainerIterator &n, INT i, MOType _mo)
{	operator [] (n)->setMO(i, _mo);	}

inline
void	extMOsSet::setClosedMO(const ContainerIterator &n, INT i, MOType _mo)
{	operator [] (n)->setMO(i+open, _mo);	}

inline
MOType	extMOsSet::getOpenMO(const ContainerIterator &n, INT i) const
{	return ((extMOsSet *) this)->operator [] (n)->getMO(i);	}

inline
MOType	extMOsSet::getClosedMO(const ContainerIterator &n, INT i) const
{	return ((extMOsSet *) this)->operator [] (n)->getMO(i+open);	}



#endif
