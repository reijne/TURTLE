//***********************************************************************
//
//	Name:			extEntry.h
//
//	Description:	
//
//	Author:			Michael Hanrath
//
//	Version:		0.0
//
//	Date:			19.03.1997
//
//
//
//
//
//***********************************************************************

#ifndef __extEntry_h
#define __extEntry_h

#include "../../../../config.h"


#include <stdarg.h>
#include <string>


#include "../../Configuration/Configuration.h"
#include "../EnergyType.h"
#include "../../MO/MOType.h"
#include "../RootEnergies.h"
#include "../Container/Tree.h"
//#include "extMOsBase.h"

class extMOsBase;
using namespace std;

class extEntry :
	public RootEnergies,
	public Tree<extMOsBase> {
public:
	extEntry(const extMOsBase *);
	extEntry(istream &s);
	extEntry(const extMOsBase *,
		INT roots, const EnergyType *e,
		INT nCSFs, const PTCIType *c,
		INT n, MOType mo ...);
	extEntry(const extMOsBase *,
		INT roots, const EnergyType *e,
		INT nCSFs, const PTCIType *c,
		INT n, MOType *mos);
	extEntry(const extMOsBase *, const RootEnergies & re, INT n, MOType *mos);
	extEntry(const extMOsBase *,
		INT	roots,
		const EnergyType *energy,
		INT nCSFs, const PTCIType *c,
		const Configuration<MOType> & conf,
		INT SelectionFlags);
	extEntry(const extMOsBase *,
		const RootEnergies &,
		const Configuration<MOType> & conf);

	// "virtual" constructor:
//	virtual extEntry * new_Set();

	void	copy(extEntry *, INT deep = 0);
	virtual extEntry *	clone(INT deep = 0);	
	
	virtual ~extEntry();

	// copy-constructor
	extEntry(const extEntry &);
	extEntry & operator = (const extEntry &);

	//------------------------------------------------------------------

	INT	getNumberOfMOs() const;

	MOType	getMO(INT i) const;
	void	setMO(INT i, MOType mo);
	const MOType *	getMOP() const;

	void	setEnergy(const RootEnergies & re);

	INT getSelectionFlags() const;
	void	setSelectionFlags(INT flag = 1);
	
	//------------------------------------------------------------------

	Configuration<MOType>	getConfiguration() const;

	//------------------------------------------------------------------

	INT	getNumberOfLeaves() const;
	void	cutTree() {}

	//--------------------------------------------------------------------------

	INT operator == (const extEntry &) const;
	INT operator != (const extEntry &) const;
	INT operator <= (const extEntry &) const;

	//----------------------------------------------------------------
	
	extEntry &	operator |= (extEntry & b);
	extEntry &	operator &= (extEntry & b);
	extEntry &	operator -= (extEntry & b);

	//------------------------------------------------------------------

	void	writeToStream(ostream & s) const;

	friend ostream& operator<<(ostream & s, const extEntry &);
	friend istream& operator>>(istream & s, extEntry &);
	
	//---------------------------------------------------------------------

private:
INT	n;						// number of MO indices
MOType	*mo;				// MO indices
INT	SelectionFlags;				// flag if configuration is selected due to near 
							// degeneration of energies in perturbation theory
};

inline
extEntry::extEntry(const extMOsBase *parent) :
		Tree<extMOsBase>(parent)
{
        n = 0;
        SelectionFlags = 0;
        mo = NULL;
}


inline
extEntry::extEntry(const extMOsBase *parent, 
	INT roots, const EnergyType *e,
	INT nCSFs, const PTCIType *c,
	INT _n, MOType _mo ...) :
		RootEnergies(roots, e, nCSFs, c), Tree<extMOsBase>(parent)
{
	n = _n;
	mo = NULL;
	SelectionFlags = 0;
	
	if ( !n )
		return;
		
	mo = new MOType[n];

va_list	ap;

	va_start(ap, _mo);
	mo[0] = _mo;
	for ( INT i=1 ; i<n ; n++ )
		mo[i] = va_arg(ap, MOType);
}


inline
extEntry::extEntry(const extMOsBase *parent, 
	INT roots, const EnergyType *e,
	INT nCSFs, const PTCIType *c,
	INT _n, MOType *mos) :
	RootEnergies(roots, e, nCSFs, c), Tree<extMOsBase>(parent)
{
	n = _n;
	mo = NULL;
	SelectionFlags = 0;
	
	if ( !n )
		return;
		
	mo = new MOType[n];
	std::memcpy(mo, mos, n*sizeof(MOType));

}

inline
extEntry::extEntry(const extMOsBase *parent,
	const RootEnergies & re, INT _n, MOType *mos) :
	RootEnergies(re), Tree<extMOsBase>(parent)
{
	n = _n;
	mo = NULL;
	SelectionFlags = 0;
	if ( !n )
		return;
		
	mo = new MOType[n];
	std::memcpy(mo, mos, n*sizeof(MOType));

}

inline
extEntry::extEntry(const extMOsBase *parent, 
	INT	roots,
	const EnergyType *e,
	INT nCSFs, const PTCIType *c,
	const Configuration<MOType> & conf,
	INT _SelectionFlags) :
	RootEnergies(roots, e, nCSFs, c), Tree<extMOsBase>(parent)
{
	n = conf.getNumberOfOpenShells() + conf.getNumberOfClosedShells();
	SelectionFlags = _SelectionFlags;
	mo = new MOType[n];
	memcpy(mo, conf.getOpenShellP(), conf.getNumberOfOpenShells()*sizeof(MOType));
	memcpy(mo+conf.getNumberOfOpenShells(),
		conf.getClosedShellP(), conf.getNumberOfClosedShells()*sizeof(MOType));
}

inline
extEntry::extEntry(const extMOsBase *parent, 
	const RootEnergies & r,
	const Configuration<MOType> & conf) :
	RootEnergies(r), Tree<extMOsBase>(parent)
{
	n = conf.getNumberOfOpenShells() + conf.getNumberOfClosedShells();
	SelectionFlags = 0;
	mo = new MOType[n];
	memcpy(mo, conf.getOpenShellP(), conf.getNumberOfOpenShells()*sizeof(MOType));
	memcpy(mo+conf.getNumberOfOpenShells(),
		conf.getClosedShellP(), conf.getNumberOfClosedShells()*sizeof(MOType));
}



inline
INT	extEntry::getNumberOfMOs() const
{	return n;	}

inline
MOType	extEntry::getMO(INT i) const
{	return mo[i];	}

inline
void	extEntry::setMO(INT i, MOType _mo)
{	mo[i] = _mo;	}


inline
INT	extEntry::getSelectionFlags() const
{	return SelectionFlags;	}

inline
void	extEntry::setSelectionFlags(INT flag)
{	SelectionFlags = flag;	}


inline
const MOType  *	extEntry::getMOP() const
{	return mo;	}


inline
void	extEntry::setEnergy(const RootEnergies & re)
{
	(RootEnergies &) *this = re;
}

inline
extEntry::extEntry(const extEntry & entry) :
	RootEnergies(entry)
{
	n = entry.getNumberOfMOs();
	SelectionFlags = entry.getSelectionFlags();
	
	mo = new MOType[n];
	memcpy(mo, entry.getMOP(), n*sizeof(MOType));
}

inline
extEntry &	extEntry::operator = (const extEntry & entry)
{
	(RootEnergies &) *this = entry;

	if ( mo )
		delete mo;
		
	n = entry.getNumberOfMOs();
	SelectionFlags = entry.getSelectionFlags();
	
	mo = new MOType[n];
	memcpy(mo, entry.getMOP(), n*sizeof(MOType));
	return *this;
}

/*
inline 
extEntry *	extEntry::
	new_Set()
{	return new extEntry();	}
*/

inline 
void	extEntry::copy(extEntry *c, INT deep)
{
	if ( deep )
	{
	// !!!!!!! deep copy missing
		*this = *c;
	}
	else
		*this = *c;
}


inline 
extEntry *
	extEntry::clone(INT deep)
{	
	extEntry * c = new extEntry(getParent());
	c->copy(this, deep);
	return c;
}


inline
INT	extEntry::getNumberOfLeaves() const
{	return 1;	}


//--------------------------------------------------------------------------


inline
INT extEntry::operator == (const extEntry &a) const
{
	if ( n!=a.getNumberOfMOs() )
		return 0;
	for ( INT i=0 ; i<n ; i++ )
		if ( mo[i]!=a.getMO(i) )
			return 0;
	return 1;
}

inline
INT extEntry::operator != (const extEntry &a) const
{
	if ( n!=a.getNumberOfMOs() )
		return 1;
	for ( INT i=0 ; i<n ; i++ )
		if ( mo[i]!=a.getMO(i) )
			return 1;
	return 0;
}

inline
INT extEntry::operator <= (const extEntry &a) const
{
	if ( n<a.getNumberOfMOs() )
		return 1;
	if ( n>a.getNumberOfMOs() )
		return 0;
	for ( INT i=0 ; i<n ; i++ )
	{
		if ( mo[i]>a.getMO(i) )
			return 0;
		if ( mo[i]<a.getMO(i) )
			return 1;
	}
	return 1;
}


#endif
