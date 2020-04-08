//***********************************************************************
//
//	Name:			extMOsDiag.h
//
//	Description:	
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

#ifndef __extMOsDiag_H
#define __extMOsDiag_H

#include "../../../../config.h"

#include "../Base/extMOsBase.h"

#include "../Base/extEntry.h"

#include "../../Configuration/Configuration.h"
#include "../../Configuration/ConfigurationSAFNr.h"

#include "TupelStructureDiag.h"
#include "../MRTreeIterator.h"

#include <iostream>

class	TupelStructureDiag;
class	extMOsSel;
class	extMOsSet;

class extMOsDiag  : 
	public extMOsBase {
public:
	extMOsDiag(istream &s);
	extMOsDiag(INT number, INT Syms, INT open, 
		INT closed, TupelStructureDiag *parent);
	~extMOsDiag();

	extMOsDiag(const extMOsDiag &);
	extMOsDiag(extMOsSel &);
	extMOsDiag(extMOsSet &, TupelStructureDiag *);


	//----------------------------------------------------------------

	const MOType *	operator [] (INT n) const;
	MOType *	operator [] (INT n);

	extEntry *& operator [] (ContainerIterator);
	extEntry * const & operator [] (ContainerIterator) const;


	const MOType *	getMOP() const;

	ConfigurationSAFNr<MOType>	getConfigurationSAFNr(INT n) const;
	
	ContainerIterator	first() const;
	void	next(ContainerIterator &) const;
	INT	isLast(ContainerIterator) const;
	
	ConfigurationSAFNr<MOType>	
		getConfigurationSAFNr(const ContainerIterator &) const;



	TupelStructureDiag * getParent() const;
		

	INT	getNumberOfElements() const;
	INT calcNumberOfLeaves() const;
	INT getNumberOfLeaves() const;
	INT	getNumberOfConfs() const;
	
	//----------------------------------------------------------------

	void	setEnergy(INT n, EnergyType e);
	void	setEnergy(EnergyType e);
	EnergyType	getEnergy(INT n = 0) const;

	//----------------------------------------------------------------

	void	setOpenMO(INT n, INT i, MOType mo, INT order=0);
	void	setClosedMO(INT n, INT i, MOType mo, INT order=0);
	
	MOType	getOpenMO(INT n, INT i, INT order=0) const;
	MOType	getClosedMO(INT n, INT i, INT order=0) const;

	INT	findOpenMO(INT order, MOType MO) const;	
	INT	findClosedMO(INT order, MOType MO) const;
	
	
	INT	getIndex(INT number, INT order) const;

	INT getSAFStart() const;
	INT getSAFStart(INT number) const;
	INT getSAFStart(INT number, INT order) const;

	//----------------------------------------------------------------

	INT	init(INT, const IrRep *);
	
	
	void	test();

//--------------------------------------------------------------------------

	void	writeToStream(ostream & s) const;

	friend ostream& operator<<(ostream & s, const extMOsDiag &);
	friend istream& operator>>(istream & s, extMOsDiag &);

//--------------------------------------------------------------------------

private:
	void	constructFromExtMOsSel(extMOsSel &);
	void	constructFromExtMOsSet(extMOsSet &);

INT	number;				// number of ext. configurations
MOType	**mo;			// array of pointers to MO lists
INT	roots;				// number of roots per configuration
INT	CSFs;				// number of CSFs per configuration
EnergyType	*energy;	// list of corresponding energies
PTCIType	*ci;		// list of corresponding ci coefficients
char	*SelectionFlags;		// array of flags if configuration is selected due to
						// near degeneration of energies in perturbation theory

INT	**MOInd;			// array of pointers to start address of 
						// specific MO (for faster access)
INT	**IndTrn;			// array of pointers to index translation table
						// 0th array contains inverse mapping

};


inline
const MOType *	extMOsDiag::getMOP() const
{	return mo[0];	}


inline
const MOType * extMOsDiag::operator [] (INT n) const
{	return mo[0]+n*total;	}

inline
MOType * extMOsDiag::operator [] (INT n)
{	return mo[0]+n*total;	}

inline
ContainerIterator	extMOsDiag::first() const
{INT nzero=0;	return ContainerIterator(nzero);	}

inline
void	extMOsDiag::next(ContainerIterator &iter) const
{	iter.i++;	}

inline
INT	extMOsDiag::isLast(ContainerIterator iter) const
{	return (iter.i>=number);	}

	
inline
ConfigurationSAFNr<MOType>	extMOsDiag::
	getConfigurationSAFNr(const ContainerIterator & iter) const
{	return getConfigurationSAFNr(iter.i);	}


inline
extEntry *& extMOsDiag::operator [] (ContainerIterator iter)
{	cout << "Error: extMOsDiag::operator [] called with " << iter.pix << endl;
	exit(1);
//	return *(operator [] (iter.i));	
}

inline
extEntry * const & extMOsDiag::operator [] (ContainerIterator iter) const
{	cout << "Error: extMOsDiag::operator [] called with " << iter.pix << endl;
	exit(1);
// 	return *(operator [] (iter.i));
}


//----------------------------------------------------------------

inline
INT	extMOsDiag::getNumberOfElements() const
{	return number;	}

inline
INT	extMOsDiag::getNumberOfLeaves() const
{	return number;  }

inline
INT	extMOsDiag::calcNumberOfLeaves() const
{	return number;  }

inline
INT	extMOsDiag::getNumberOfConfs() const
{	return number;	}

//----------------------------------------------------------------

inline
void	extMOsDiag::setEnergy(INT n, EnergyType e)
{	energy[n] = e;	}

inline
void	extMOsDiag::setEnergy(EnergyType e)
{	energy[0] = e;	}

inline
EnergyType	extMOsDiag::getEnergy(INT n) const
{	return energy[n];	}

//----------------------------------------------------------------

inline
void	extMOsDiag::setOpenMO(INT n, INT i, MOType _mo, INT order)
{	*(mo[order]+n*total+i) = _mo;	}

inline
void	extMOsDiag::setClosedMO(INT n, INT i, MOType _mo, INT order)
{	*(mo[order]+n*total+open+i) = _mo;	}

inline
MOType	extMOsDiag::getOpenMO(INT n, INT i, INT order) const
{	return *(mo[order]+n*total+i);	}

inline
MOType	extMOsDiag::getClosedMO(INT n, INT i, INT order) const
{	return *(mo[order]+n*total+open+i);	}

inline
INT	extMOsDiag::findOpenMO(INT order, MOType MO) const
{
	if ( MO>maxMO )
		return -1;
	else
		return MOInd[order][MO];	
}

inline
INT	extMOsDiag::findClosedMO(INT order, MOType MO) const
{	return findOpenMO(order+open, MO);	}	

inline
INT	extMOsDiag::getIndex(INT number, INT order) const
{
	if ( !order )
		return number;
	else
		return	IndTrn[order][number];
}


inline
INT extMOsDiag::getSAFStart() const
{	return extMOsBase::getSAFStart();	}

inline
INT extMOsDiag::getSAFStart(INT number) const
{	return extMOsBase::getSAFStart(number);	}

inline
INT	extMOsDiag::getSAFStart(INT number, INT order) const
{	return SAFStart + getIndex(number, order)*SAFInc;	}


#endif
