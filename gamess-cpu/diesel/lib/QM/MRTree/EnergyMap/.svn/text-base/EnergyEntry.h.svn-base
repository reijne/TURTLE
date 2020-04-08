//***********************************************************************
//
//	Name:			EnergyEntry.h
//
//	Description:	
//
//	Author:			Michael Hanrath
//
//	Version:		0.0
//
//	Date:			19.04.1997
//
//
//
//
//***********************************************************************

#ifndef __EnergyEntry_h
#define __EnergyEntry_h

#include "../../../../config.h"

#include "../EnergyType.h"
#include "../PTCIType.h"

class EnergyEntry {
public:
	EnergyEntry(INT roots, INT safs);
	~EnergyEntry();

enum SelectionMode	{ AmountSquareOfSum, SumOfAmountSquares };
enum Status	{ NotSelected = 0, Selected = 1, SelectionFlagserated = 2, Random = 3 };
static const EnergyType	DegenerationThreshold;


	INT	getNumberOfRoots() const;
	INT getNumberOfSAFs() const;
	EnergyType *getNominatorP() const;
	EnergyType *getDenominatorP() const;
	
	EnergyType	getEnergy(INT root) const;
	const EnergyType	*getEnergyP() const;

	PTCIType	getCICoef(INT root, INT CSF) const;
	const PTCIType	*getCICoefP() const;

	Status	select(
		EnergyType threshold,
		SelectionMode mode);


private:
INT	nRoots;					// number of roots
INT	nSAFs;					// number of SAFs
EnergyType	*nominator;		// pointer to nominators of energy expression
EnergyType	*denominator;	// pointer to denominators of energy expression
EnergyType	*dE;			// pointer to PT energies
PTCIType	*ci;			// pointer to estimated ci vectors
};


#include <string>

inline
EnergyEntry::EnergyEntry(INT roots, INT safs)
{
	nSAFs = safs;
	nRoots = roots;
	nominator = new EnergyType[roots*safs];
	denominator = new EnergyType[roots*safs];
	std::memset(nominator, 0, roots*safs*sizeof(EnergyType));
	dE = NULL;
	ci = NULL;
}


inline
EnergyEntry::~EnergyEntry()
{
	if ( nominator )
		delete[] nominator;
	if ( denominator)
		delete[] denominator;
	if ( dE )
		delete[] dE;
	if ( ci )
		delete[] ci;
}


inline
INT	EnergyEntry::getNumberOfRoots() const
{	return	nRoots;	}

inline
INT	EnergyEntry::getNumberOfSAFs() const
{	return	nSAFs;	}

inline
EnergyType	*EnergyEntry::getNominatorP() const
{	return nominator;	};

inline
EnergyType	*EnergyEntry::getDenominatorP() const
{	return denominator;	};

inline
EnergyType	EnergyEntry::getEnergy(INT root) const
{	return	dE[root];	}

inline
const EnergyType	*EnergyEntry::getEnergyP() const
{	return	dE;	}

inline
PTCIType	EnergyEntry::getCICoef(INT root, INT CSF) const
{	return	ci[root*nSAFs + CSF];	}

inline
const PTCIType	*EnergyEntry::getCICoefP() const
{	return	ci;	}




#endif
