//***********************************************************************
//
//	Name:	Configuration.h
//
//	Description:	implements configuration handling
//					based on second quantization
//
//	Author:			Michael Hanrath
//
//	Version:		0.0
//
//	Date:			21.06.1996
//
//
//
//
//
//***********************************************************************

#ifndef __CONFIGURATION_H
#define __CONFIGURATION_H

#include "../../../config.h"

#include <iostream>
#include <stdarg.h>

#include "../../Math/etc/MathObject.h"

#include "ConfigurationGlobals.h"
#include "../MO/MOType.h"
#include "../MO/GeneralizedMO.h"
#include "../MO/Iterators/MOListIterator.h"

class Excitation;
class	TupelStructureSel;
class	MRMOs;

template <class TMOType>
class Configuration;

template <class TMOType>
istream& operator>> (istream & s, Configuration<TMOType> &);

template <class TMOType>
class Configuration : public MathObject {
public:
//friend class Configuration<MOType>;
//friend class Configuration<GeneralizedMO>;
	Configuration();
	Configuration(istream &s);
	Configuration(INT NumberOfOpenShells, INT NumberOfClosedShells);
	Configuration(INT NumberOfOpenShells, INT NumberOfClosedShells, 
		const TMOType *pOpen, const TMOType *pClosed);
	Configuration(INT NumberOfOpenShells, INT NumberOfClosedShells, 
		const TMOType *pShells);
	Configuration(INT NumberOfOpenShells, INT NumberOfClosedShells, 
		INT mo1);
	Configuration(INT NumberOfOpenShells, INT NumberOfClosedShells, 
		GeneralizedMO mo1);
	Configuration(INT NumberOfOpenShells, INT NumberOfClosedShells, 
		INT mo1, INT mo2);
	Configuration(INT NumberOfOpenShells, INT NumberOfClosedShells, 
		GeneralizedMO mo1, GeneralizedMO mo2);
/*
	varargs do not work on SGI with GNU C/C++ 2.8.0
	Configuration(INT NumberOfOpenShells, INT NumberOfClosedShells, 
		INT mo, ...);
	Configuration(INT NumberOfOpenShells, INT NumberOfClosedShells, 
		GeneralizedMO mo, ...);
*/	
	
	// "virtual" constructor:
/*	virtual Configuration<TMOType> * new_Set();
	
	void	copy(Configuration<TMOType> *, INT deep = 0);
	virtual Configuration<TMOType> *	clone(INT deep = 0);
*/
	// Copy-Konstruktor
	Configuration(const Configuration<TMOType> &);
	
	// Zuweisungs-Operator
	Configuration & operator = (const Configuration<TMOType> &);

	// type conversion
	operator Configuration<GeneralizedMO> ();
	Configuration(const Configuration<GeneralizedMO> &,
		const MOListIterator &);

	Configuration<TMOType> & operator * () { return *this; }
//--------------------------------------------------------------------------

	INT	annihilate(TMOType mo);
	INT	create(TMOType mo);
	INT	create(Configuration<TMOType> conf);
	void	create(const MOListIterator & molist);
	INT	excite(TMOType from, TMOType to);
	INT	excite(const Excitation &);
	
	INT	deleteSingleMO(TMOType MO);
	INT	deleteDoubleMO(TMOType MO);

//--------------------------------------------------------------------------

	static void	calcInteraction(
		const Configuration<TMOType> & a,
		const Configuration<TMOType> & b, 
		Configuration<TMOType> & da, Configuration<TMOType> & db,
		INT *iposa, INT *iposb);

	// general method
	static INT calcExcitationOrder(
		const Configuration<TMOType> & a,
		const Configuration<TMOType> & b);

	// much faster method, requires a to have more or equal electrons than b
	// returns maxOrder+1 for any excitation order greater maxOrder
	static INT calcExcitationOrderFast(
		const Configuration<TMOType> & a,
		const Configuration<TMOType> & b, INT maxOrder = 2);


	
//--------------------------------------------------------------------------

	void	split(const MRMOs *,
		Configuration<MOType> &internal, Configuration<MOType> &external) const;
	
	INT	getNumberOfExternals(const MRMOs *) const;

//--------------------------------------------------------------------------

	//	only for MOType
	IrRep	calcIrRep(const MOIrReps & moirreps) const;

	//	only for GeneralizedMO
	void	setIrRep_Space();
	INT	getNthExtPosOpen(INT) const;

//--------------------------------------------------------------------------

	const TMOType	*getOpenShellP() const;
	const TMOType	*getClosedShellP() const;

	const TMOType &	getOpenShell(INT n) const;
	const TMOType &	getClosedShell(INT n) const;
	const TMOType &	getShell(INT n) const;

	void	setOpenShell(INT i, TMOType mo);
	void	setClosedShell(INT i, TMOType mo);

	TMOType &	getOpenShellRef(INT n);
	TMOType &	getClosedShellRef(INT n);

	INT	getNumberOfOpenShells() const;
	INT	getNumberOfClosedShells() const;
	INT	getNumberOfElectrons() const;
	
	void	insertOpenMO(TMOType mo);
	void	insertClosedMO(TMOType mo);

	void	appendOpenShell(TMOType mo);
	void	appendClosedShell(TMOType mo);

	void	append(Configuration<TMOType> conf);
	
//--------------------------------------------------------------------------
	
	// faster than sort (supposes two parts (internal and external)
	// are already sorted among themselves)
	void	sortMerge(INT openExtStart, INT closedExtStart);
	void	sortMerge(INT openExtStart, INT closedExtStart, INT *ipos);

	void	clear();

//--------------------------------------------------------------------------

/*	friend	INT operator == <> (const Configuration<TMOType> &,
		const Configuration<TMOType> &);
	friend	INT operator != <> (const Configuration<TMOType> &,
		const Configuration<TMOType> &);
	friend	INT operator <= <> (const Configuration<TMOType> &,
		const Configuration<TMOType> &);
*/
//--------------------------------------------------------------------------

	Configuration<TMOType> operator - 
		(const Configuration<TMOType> &);

	Configuration<TMOType> &operator &= 
		(const Configuration<TMOType> &);



	Configuration<TMOType>  operator & (const Configuration &);
	Configuration<TMOType> & operator += (const Configuration &);
	Configuration<TMOType> & operator -= (const Configuration &);


//--------------------------------------------------------------------------

	void	writeToStream(ostream & s) const;

//--------------------------------------------------------------------------

	INT	check() const;
	void	sort();

//--------------------------------------------------------------------------

          template <class TT> friend      ostream& operator<< (ostream & s, const Configuration<TT> &);
          friend  istream& operator>> <> (istream & s, Configuration<TMOType> &);

private:
	INT	insertMO(INT & shells, TMOType *p, TMOType mo);
	INT	deleteMO(INT & shells, TMOType *p, TMOType mo);
	
	
	
INT	openShells;						// number of open shells
INT	closedShells;					// number of closed shells
TMOType	pOpen[MAXELECTRONS];		// open shells
TMOType	pClosed[MAXELECTRONS];		// closed shells

};

template <class TMOType>
inline
INT	operator != (const Configuration<TMOType> & conf1, const Configuration<TMOType> & conf2)
{	return !(conf1==conf2);	}


INT	operator == (const Configuration<MOType> &, const Configuration<MOType> &);
INT	operator <= (const Configuration<MOType> &, const Configuration<MOType> &);
INT	operator < (const Configuration<MOType> &, const Configuration<MOType> &);
/*
template <class TMOType>	INT operator == (const Configuration<TMOType> &,
		const Configuration<TMOType> &);
template <class TMOType>	INT operator != (const Configuration<TMOType> &,
		const Configuration<TMOType> &);
template <class TMOType>	INT operator <= (const Configuration<TMOType> &,
		const Configuration<TMOType> &);
*/


//istream& operator>>(istream & s, Configuration<MOType> &conf);

//=========================================================================

#include <string>

template <class TMOType>
inline
Configuration<TMOType>::Configuration(const Configuration<TMOType> & conf)
{
	openShells = conf.getNumberOfOpenShells();
	closedShells = conf.getNumberOfClosedShells();
	memcpy(pOpen, conf.getOpenShellP(), openShells*sizeof(TMOType));
	memcpy(pClosed, conf.getClosedShellP(), closedShells*sizeof(TMOType));
}
	
template <class TMOType>
inline
Configuration<TMOType> & Configuration<TMOType>::operator = 
	(const Configuration<TMOType> & conf)
{
	openShells = conf.getNumberOfOpenShells();
	closedShells = conf.getNumberOfClosedShells();
	memcpy(pOpen, conf.getOpenShellP(), openShells*sizeof(TMOType));
	memcpy(pClosed, conf.getClosedShellP(), closedShells*sizeof(TMOType));
	return *this;
}



template <class TMOType>
inline
Configuration<TMOType>::Configuration()
{	closedShells = openShells = 0;	}

template <class TMOType>
inline
Configuration<TMOType>::Configuration(INT NumberOfOpenShells, INT NumberOfClosedShells)
{
	openShells = NumberOfOpenShells;
	closedShells = NumberOfClosedShells;	
}


template <class TMOType>
inline
Configuration<TMOType>::Configuration(INT NumberOfOpenShells, INT NumberOfClosedShells,
	const TMOType *_pOpen, const TMOType *_pClosed)
{
	openShells = NumberOfOpenShells;
	closedShells = NumberOfClosedShells;
	memcpy(pOpen, _pOpen, openShells*sizeof(TMOType));
	memcpy(pClosed, _pClosed, closedShells*sizeof(TMOType));
}


template <class TMOType>
inline
Configuration<TMOType>::Configuration(INT NumberOfOpenShells, INT NumberOfClosedShells,
	const TMOType *_pShells)
{
	openShells = NumberOfOpenShells;
	closedShells = NumberOfClosedShells;
	memcpy(pOpen, _pShells, openShells*sizeof(TMOType));
	memcpy(pClosed, _pShells + openShells, 
		closedShells*sizeof(TMOType));
}


/*
template <class TMOType>
inline 
Configuration<TMOType> *	Configuration<TMOType>::
	new_Set()
{	return new Configuration<TMOType>();	}



template <class TMOType>
inline 
void	Configuration<TMOType>::copy(Configuration<TMOType> *c, INT deep)
{
	if ( deep )
	{
	// !!! deep copy to implement
		*this = *c;
	}
	else
		*this = *c;
}


template <class TMOType>
inline 
Configuration<TMOType> *	Configuration<TMOType>::clone(INT deep)
{	
	Configuration<TMOType> * c = new Configuration<TMOType>();
	c->copy(this, deep);
	return c;
}
*/



//--------------------------------------------------------------------------


template <class TMOType>
inline
const TMOType	*Configuration<TMOType>::getOpenShellP() const
{	return	pOpen;	}

template <class TMOType>
inline
const TMOType	*Configuration<TMOType>::getClosedShellP() const
{	return	pClosed;	}


template <class TMOType>
inline
const TMOType &	Configuration<TMOType>::getOpenShell(INT n) const
{	return	pOpen[n];	}

template <class TMOType>
inline
const TMOType &	Configuration<TMOType>::getClosedShell(INT n) const
{	return	pClosed[n];	}

template <class TMOType>
inline
const TMOType &	Configuration<TMOType>::getShell(INT n) const
{	if ( n<openShells )
		return	pOpen[n];
	else
		return	pClosed[n-openShells];
}


template <class TMOType>
inline
TMOType &	Configuration<TMOType>::getOpenShellRef(INT n)
{	return	pOpen[n];	}

template <class TMOType>
inline
TMOType &	Configuration<TMOType>::getClosedShellRef(INT n)
{	return	pClosed[n];	}


template <class TMOType>
inline
INT	Configuration<TMOType>::getNumberOfOpenShells() const
{	return	openShells;	}

template <class TMOType>
inline
INT	Configuration<TMOType>::getNumberOfClosedShells() const
{	return	closedShells;	}

template <class TMOType>
inline
INT	Configuration<TMOType>::getNumberOfElectrons() const
{	return	openShells + (closedShells << 1);	}

template <class TMOType>
inline
void	Configuration<TMOType>::appendOpenShell(TMOType mo)
{	pOpen[openShells++] = mo;	}

template <class TMOType>
inline
INT	Configuration<TMOType>::create(Configuration<TMOType> conf)
{
	for ( INT i=0 ; i<conf.getNumberOfOpenShells() ; i++ )
		if ( create(conf.getOpenShell(i)) )
			return 1;

	for ( INT i=0 ; i<conf.getNumberOfClosedShells() ; i++ )
		insertClosedMO(conf.getClosedShell(i));
	return 0;
}

template <class TMOType>
inline
void	Configuration<TMOType>::append(Configuration<TMOType> conf)
{
	for ( INT i=0 ; i<conf.getNumberOfOpenShells() ; i++ )
		pOpen[openShells++] = conf.getOpenShell(i);

	for ( INT i=0 ; i<conf.getNumberOfClosedShells() ; i++ )
		pClosed[closedShells++] = conf.getClosedShell(i);
}

template <class TMOType>
inline
void	Configuration<TMOType>::appendClosedShell(TMOType mo)
{	pClosed[closedShells++] = mo;	}

template <class TMOType>
inline
void	Configuration<TMOType>::clear()
{	closedShells = openShells = 0;	}

template <class TMOType>
inline
void	Configuration<TMOType>::setOpenShell(INT i, TMOType mo)
{	pOpen[i] = mo;	}

template <class TMOType>
inline
void	Configuration<TMOType>::setClosedShell(INT i, TMOType mo)
{	pClosed[i] = mo;	}




#endif


