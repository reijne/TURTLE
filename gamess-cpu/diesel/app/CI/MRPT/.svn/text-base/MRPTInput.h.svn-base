//***********************************************************************
//
//	Name:			MRPTInput.h
//
//	Description:	stores configuration input for 
//					individually selecting MR-CI
//
//	Author:			Michael Hanrath
//
//	Version:		0.0
//
//	Date:			05.02.1997
//
//
//
//
//
//***********************************************************************

#ifndef __MRPTInput_H
#define __MRPTInput_H

#include <iostream>

#include "../../../config.h"

#include "../../../lib/Container/String.h"

#include "../../../lib/QM/RepresentationMatrices/DataTypes.h"

#include "../../../lib/QM/MRTree/EnergyType.h"
#include "../../../lib/QM/Configuration/Configuration.h"
#include "../../../lib/QM/IO/Fortran/Fort31File.h"
#include "../../../lib/QM/Configuration/ConfigurationSet.h"
#include "../Selector/EnergyMap.h"

#include "../../../lib/QM/MO/MOMapping.h"

#include "../../../lib/QM/MRCIMatrix/MRMPH0Matrix.h"

class MRMOs;
class MORestriction;
class MOEquivalence;
class MOStatistics;

class MRPTInput {
public:
	MRPTInput();
	~MRPTInput();

friend class yyFlexLexer;


	// copy-constructor
	MRPTInput(const MRPTInput &);

	MRPTInput & operator = (const MRPTInput &);

	//--------------------------------------------------------------------------

	Fort31File	getMOIntegralFile() const;

	INT	getNumberOfSelectionThresholds() const;
	EnergyType	getSelectionThreshold(INT i) const;
	const char	*getSelectionThresholdString(INT i) const;
	const EnergyType	*getSelectionThresholdP() const;


	const MRMOs	*getMRMOs() const;


	INT	getNumberOfRoots() const;
	INT getRootNumber(INT i) const;
	const INT * getRootNumberP() const;

	RootType getRootEnergy(INT i) const;
	const RootType* getRootEnergyP() const;

	void	setVerbosity(const char *);


	INT	getCalcMP3() const;


	MRMPH0Matrix<double, double>::ProjectionMode getProjectionMode() const;
	String	getInhomogenityThreshold() const;
	
	
	//--------------------------------------------------------------------------

	friend ostream& operator<<(ostream & s, const MRPTInput &);
	friend istream& operator>>(istream & s, MRPTInput &);

	//--------------------------------------------------------------------------

protected:
Fort31File	MOIntegralFile;			// MO integral filename
INT	NumberOfElectrons;

EnergyType	*SelectionThreshold;	// threshold for selection procedure
char	**SelectionThresholdString;	// same, but as a character string
INT	nSelectionThresholds;			// number of thresholds

INT	NumberOfRoots;
INT	NumberOfRootsE;
INT	*rootNumbers;
RootType	*rootEnergies;
MRMOs	*mrmos;

MRMPH0Matrix<double, double>::ProjectionMode projectionMode;
String	inhomogenityThreshold;
INT	calcMP3;

private:
};



inline
Fort31File	MRPTInput::getMOIntegralFile() const
{	return MOIntegralFile;	}

inline
INT	MRPTInput::getNumberOfSelectionThresholds() const
{	return nSelectionThresholds;	}

inline
EnergyType	MRPTInput::getSelectionThreshold(INT i) const
{	return SelectionThreshold[i];	}

inline
const char	*MRPTInput::getSelectionThresholdString(INT i) const
{	return SelectionThresholdString[i];	}

inline
const EnergyType	*MRPTInput::getSelectionThresholdP() const
{	return SelectionThreshold;	}



inline
const MRMOs	*MRPTInput::getMRMOs() const
{	return mrmos;	}

inline
INT	MRPTInput::getNumberOfRoots() const
{	return NumberOfRoots;	}

inline
INT	MRPTInput::getRootNumber(INT i) const
{	return rootNumbers[i];	}


inline
const INT * MRPTInput::getRootNumberP() const
{	return rootNumbers;	}


inline
RootType MRPTInput::getRootEnergy(INT i) const
{	return rootEnergies[i];	}

inline
const RootType* MRPTInput::getRootEnergyP() const
{	return rootEnergies;	}

inline
MRMPH0Matrix<double, double>::ProjectionMode MRPTInput::getProjectionMode() const
{	return projectionMode;	}

inline
String	MRPTInput::getInhomogenityThreshold() const
{	return inhomogenityThreshold;	}

inline
INT	MRPTInput::getCalcMP3() const
{	return	calcMP3;	}


#endif
