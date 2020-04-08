//***********************************************************************
//
//	Name:			DiagInput.h
//
//	Description:	stores configuration input for 
//					individually selecting MR-CI
//
//	Author:			Michael Hanrath
//
//	Version:		0.0
//
//	Date:			13.09.1997
//
//
//
//
//
//***********************************************************************


#ifndef __DiagInput_H
#define __DiagInput_H

#include <iostream>

#include "../../../config.h"

#include "../../../lib/Container/String.h"


#include "../../../lib/QM/RepresentationMatrices/CIVectors.h"
#include "../../../lib/QM/MRTree/EnergyType.h"
#include "../../../lib/QM/Configuration/Configuration.h"
#include "../../../lib/QM/IO/Fortran/Fort31File.h"
#include "../../../lib/QM/Configuration/ConfigurationSet.h"
#include "../../../lib/QM/Davidson/DavidsonCI.h"


class MRMOs;

class DiagInput {
public:
	DiagInput(INT autoDetectFort31);
	~DiagInput();

friend class yyFlexLexer;

typedef double	VectorType;


enum Precision { undefinedPrec, floatPrec, doublePrec };

	// copy-constructor
	DiagInput(const DiagInput &);

	//--------------------------------------------------------------------------

	Fort31File	getMOIntegralFile() const;

	VectorType	getRefThreshold() const;
	VectorType	getPTRefThreshold() const;
	EnergyType		getConvergenceEnergyChange() const;
	VectorType	getConvergenceEigenvectorChange() const;
        String          getConfTreeFileName() const;
	INT				getMaxIters() const;
	INT				getMaxStorageMem() const;

	Precision	getPrecision() const;

	ConfigurationSet &	getRefConfSet();
	ConfigurationSet &	getPTRefConfSet();

	INT	getNumberOfRoots() const;
	INT getRootNumber(INT i) const;
	const INT *getRootNumbersP() const;

	INT	getRootHoming() const;
	INT	getStorePTEnergy() const;
	INT	getStorePTCoef() const;

	DavidsonCI<double, double>::IterationMode	getIterationMode() const;
	
	void	setVerbosity(const char *);
	
	//--------------------------------------------------------------------------

	friend ostream& operator<<(ostream & s, const DiagInput &);
	friend istream& operator>>(istream & s, DiagInput &);

	//--------------------------------------------------------------------------

protected:
Fort31File	MOIntegralFile;			// MO integral filename
									

VectorType	RefThreshold;				// threshold for reference selection
VectorType	PTRefThreshold;				// threshold for 0th order wave function

EnergyType	ConvergenceEnergyChange;		// convergence threshold (energy)

VectorType	ConvergenceEigenvectorChange;	// convergence threshold (eigen vector)

String          ConfTreeFileName;                     // filename for ConfTree.dat file


INT	rootHoming;								// flag if root homing is active

INT	storePTEnergy;							// flag if PT energies are stored in
											// configuration tree

INT	storePTCoef;							// flag if PT coefs are stored in
											// configuration tree

INT MaxIters;								// maximum number of iterations
//LONG_INT	MaxStorageMem;					// maximum memory consumption 
INT      MaxStorageMem;                                  // maximum memory consumption
											// in conventional (non direct) mode

DavidsonCI<double, double>::IterationMode	iterationMode;

INT	NumberOfRoots;
INT	*rootNumbers;

Precision	precision;


static const VectorType refThreshDefault;
static const VectorType PTrefThreshDefault;
static const EnergyType ConvergenceEnergyChangeDefault;
static const VectorType ConvergenceEigenvectorChangeDefault;
static const String ConfTreeFileNameDefault;
static const INT MaxItersDefault;
static const DavidsonCI<double, double>::IterationMode	iterationModeDefault;
private:
};



inline
Fort31File	DiagInput::getMOIntegralFile() const
{	return MOIntegralFile;	}

inline
DiagInput::VectorType	DiagInput::getRefThreshold() const
{	return RefThreshold;	}

inline
DiagInput::VectorType	DiagInput::getPTRefThreshold() const
{	return PTRefThreshold;	}

inline
EnergyType	DiagInput::getConvergenceEnergyChange() const
{	return ConvergenceEnergyChange;	}

inline
DiagInput::VectorType	DiagInput::getConvergenceEigenvectorChange() const
{	return ConvergenceEigenvectorChange;	}

inline 
String DiagInput::getConfTreeFileName() const
{       return ConfTreeFileName; }

inline
INT	DiagInput::getMaxIters() const
{	return MaxIters;	}

inline
INT	DiagInput::getMaxStorageMem() const
{	return MaxStorageMem;	}

inline
INT	DiagInput::getNumberOfRoots() const
{	return NumberOfRoots;	}

inline
INT	DiagInput::getRootNumber(INT i) const
{	return rootNumbers[i];	}

inline
const INT	*DiagInput::getRootNumbersP() const
{	return rootNumbers;	}

inline
INT	DiagInput::getRootHoming() const
{	return	rootHoming;	}

inline
INT	DiagInput::getStorePTEnergy() const
{	return	storePTEnergy;	}

inline
INT	DiagInput::getStorePTCoef() const
{	return	storePTCoef;	}

inline
DavidsonCI<double, double>::IterationMode	DiagInput::getIterationMode() const
{	return	iterationMode;	}

inline
DiagInput::Precision	DiagInput::getPrecision() const
{	return	precision;	}

#endif
