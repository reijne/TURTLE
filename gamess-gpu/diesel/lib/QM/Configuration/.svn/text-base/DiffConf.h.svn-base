//***********************************************************************
//
//	Name:	DiffConf.h
//
//	Description:	general representation of excitations
//
//	Author:	Michael Hanrath
//
//	Version:	0.0
//
//	Date:	11.06.1996
//
//
//
//
//
//***********************************************************************


#ifndef __DIFFCONF_H
#define __DIFFCONF_H

#include "../../../config.h"

#include <iostream>
#include <stdarg.h>

#include "../../Math/etc/MathObject.h"

#include "Configuration.h"
#include "ConfigurationGlobals.h"
#include "../MO/MOType.h"
#include "../MO/GeneralizedMO.h"
#include "../MO/Iterators/MOListIterator.h"


template <class TMOType>
class DiffConf : public MathObject {
public:
//friend class DiffConf<MOType>;
//friend class DiffConf<GeneralizedMO>;
	DiffConf() {}
//	DiffConf(Configuration<TMOType> & from, Configuration<TMOType> & to);
	~DiffConf() {}
	
	
	// copy constructor
	DiffConf(const DiffConf &);
	
	// assignment operator
	DiffConf & operator = (const DiffConf &);

	// type conversion
	operator DiffConf<GeneralizedMO> ();
	DiffConf(const DiffConf<GeneralizedMO> &, const MOListIterator &);

//--------------------------------------------------------------------------


	INT	getOrder() const;

	void	calcDiffConf(const Configuration<TMOType> & a, const Configuration<TMOType> & b,
		INT append = 0);


	// (only for external MOs !)
	void	addExternal(const Configuration<TMOType> & a,
		const Configuration<TMOType> & b);

	INT	addExternalTo(INT NumberOfOpenShells, INT NumberOfClosedShells,
		const TMOType *_pShells);

//	void	addExternal(INT nFrom, *MOType moFrom, INT nTo, *MOType moTo);

//--------------------------------------------------------------------------

	void	setSame(const Configuration<TMOType> & conf);

	const Configuration<TMOType> &	getSame() const;
	const Configuration<TMOType> &	getFrom() const;
	const Configuration<TMOType> &	getTo() const;

	void	setToOpen(INT i, TMOType mo);
	void	setToClosed(INT i, TMOType mo);
	void	setFromOpen(INT i, TMOType mo);
	void	setFromClosed(INT i, TMOType mo);

	INT	getPosFrom(INT i) const;
	INT	getPosTo(INT i) const;
	
	const INT	*getPosFromP() const;
	const INT	*getPosToP() const;

	INT	getOpenShellsFrom() const;
	INT	getOpenShellsTo() const;


	// faster than sort (supposes two parts (internal and external)
	// are already sorted among themselves)
	void	updatePos(Configuration<TMOType> & conf, INT *ipos);

//--------------------------------------------------------------------------

	void	clear();
	
//--------------------------------------------------------------------------

	Configuration<MOType>	simplify();
	
//--------------------------------------------------------------------------

//--------------------------------------------------------------------------

	template <class TT> friend ostream& operator<< (ostream & s, const DiffConf<TT> &);

//--------------------------------------------------------------------------




private:



Configuration<TMOType>	same;	// identically occupied MOs in a and b
Configuration<TMOType>	from;	// missing MOs in b from a
Configuration<TMOType>	to;		// missing MOs in a from b

INT	openShellsFrom;				// number of open shells in
INT	openShellsTo;				// original configurations

INT posFrom[MAXELECTRONS];		// position of interacting open shell in
INT	posTo[MAXELECTRONS];		// original configurations a and b
								// (needed for calculation of q-case)
};


template <class TMOType>	INT operator == (const DiffConf<TMOType> &, const DiffConf<TMOType> &);
template <class TMOType>	INT operator != (const DiffConf<TMOType> &, const DiffConf<TMOType> &);


//template <class TMOType>	ostream& operator<<(ostream & s, const DiffConf<TMOType> &);



//=========================================================================


#include <string>


template <class TMOType>
inline
void	DiffConf<TMOType>::setSame(const Configuration<TMOType> & conf)
{
	same = conf;
	from.clear();
	to.clear();
	openShellsFrom = openShellsTo = same.getNumberOfOpenShells();
}

template <class TMOType>
inline
DiffConf<TMOType>::DiffConf(const DiffConf & exc)
{
	same = exc.getSame();

	from = exc.getFrom();
	memcpy(posFrom, exc.getPosFromP(), from.getNumberOfOpenShells()*sizeof(INT));
	openShellsFrom = exc.getOpenShellsFrom();
		
	to = exc.getTo();
	memcpy(posTo, exc.getPosToP(), to.getNumberOfOpenShells()*sizeof(INT));
	openShellsTo = exc.getOpenShellsTo();
}
	
template <class TMOType>
inline
DiffConf<TMOType> & DiffConf<TMOType>::operator = (const DiffConf<TMOType> & exc)
{
	same = exc.getSame();

	from = exc.getFrom();
	memcpy(posFrom, exc.getPosFromP(), from.getNumberOfOpenShells()*sizeof(INT));
	openShellsFrom = exc.getOpenShellsFrom();
		
	to = exc.getTo();
	memcpy(posTo, exc.getPosToP(), to.getNumberOfOpenShells()*sizeof(INT));
	openShellsTo = exc.getOpenShellsTo();

	return *this;
}

template <class TMOType>
inline
const Configuration<TMOType> &	DiffConf<TMOType>::getSame() const
{	return	same;	}

template <class TMOType>
inline
const Configuration<TMOType> &	DiffConf<TMOType>::getFrom() const
{	return	from;	}

template <class TMOType>
inline
const Configuration<TMOType> &	DiffConf<TMOType>::getTo() const
{	return	to;	}

template <class TMOType>
inline
INT	DiffConf<TMOType>::getPosFrom(INT i) const
{	return	posFrom[i];	}

template <class TMOType>
inline
INT	DiffConf<TMOType>::getPosTo(INT i) const
{	return	posTo[i];	}

template <class TMOType>
inline
const INT *	DiffConf<TMOType>::getPosFromP() const
{	return	posFrom;	}

template <class TMOType>
inline
const INT *	DiffConf<TMOType>::getPosToP() const
{	return	posTo;	}

template <class TMOType>
inline
INT	DiffConf<TMOType>::getOpenShellsFrom() const
{	return openShellsFrom;	}

template <class TMOType>
inline
INT	DiffConf<TMOType>::getOpenShellsTo() const
{	return openShellsTo;	}

template <class TMOType>
inline
void	DiffConf<TMOType>::clear()
{
	same.clear();
	from.clear();
	to.clear();
	openShellsFrom = openShellsTo = 0;
}

template <class TMOType>
inline
void	DiffConf<TMOType>::setToOpen(INT i, TMOType mo)
{	to.setOpenShell(i, mo);	}

template <class TMOType>
inline
void	DiffConf<TMOType>::setToClosed(INT i, TMOType mo)
{	to.setClosedShell(i, mo);	}

template <class TMOType>
inline
void	DiffConf<TMOType>::setFromOpen(INT i, TMOType mo)
{	from.setOpenShell(i, mo);	}

template <class TMOType>
inline
void	DiffConf<TMOType>::setFromClosed(INT i, TMOType mo)
{	from.setClosedShell(i, mo);	}


template <class TMOType>
inline
INT	DiffConf<TMOType>::getOrder() const
{
	return Configuration<TMOType>::calcExcitationOrder(from, to);
}


#endif
