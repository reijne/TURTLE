//***********************************************************************
//
//	Name:			RefSelInput.h
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

#ifndef __RefSelInput_H
#define __RefSelInput_H

#include <iostream>

#include "../../../config.h"

#include "../../../lib/Container/String.h"
#include "../../../lib/Container/SLList.h"

#include "../../../lib/QM/RepresentationMatrices/DataTypes.h"

#include "../../../lib/QM/MRTree/EnergyType.h"
#include "../../../lib/QM/Configuration/Configuration.h"
#include "../../../lib/QM/IO/Fortran/Fort31File.h"
#include "../../../lib/QM/Configuration/ConfigurationSet.h"
#include "../Selector/EnergyMap.h"

#include "../../../lib/QM/MO/MOMapping.h"
#include "../../../lib/QM/Symmetry/DegenPointgroups/RefSelMode.h"


class MRMOs;
class MORestriction;
class MOEquivalence;
class MOStatistics;

class RefSelInput {
public:
	RefSelInput();
	~RefSelInput();

friend class yyFlexLexer;



	// copy-constructor
	RefSelInput(const RefSelInput &);

	RefSelInput & operator = (const RefSelInput &);

	//--------------------------------------------------------------------------

	Fort31File	getMOIntegralFile() const;


	
	INT	getFirstGuessConfs() const;
	
	double	getReferenceThreshold() const;


	INT	getNumberOfSelectionThresholds() const;
	EnergyType	getSelectionThreshold(INT i) const;
	const char	*getSelectionThresholdString(INT i) const;
	const EnergyType	*getSelectionThresholdP() const;


	const MOEquivalence	*getMOEquivalence() const;
	void	setMOEquivalence(const char *, INT auto = 0);


	const MRMOs	*getMRMOs() const;


	RefSelMode::Mode getRefSelMode() const
	{	return refSelMode;	}


	void	setVerbosity(const char *);

	SLList<INT>	getIrReps()
	{	return irreps;	}

	//--------------------------------------------------------------------------

	friend ostream& operator<<(ostream & s, const RefSelInput &);
	friend istream& operator>>(istream & s, RefSelInput &);

	//--------------------------------------------------------------------------

protected:
Fort31File	MOIntegralFile;			// MO integral filename


INT	firstGuessConfs;				// number of configurations used in
									// first guess minus #roots


INT	autoEquiv;						// determine degenerated MOs
									// automatically


RefSelMode::Mode	refSelMode;

double	ReferenceThreshold;

MRMOs	*mrmos;

MOEquivalence	*moEquivalence;
SLList<INT>	irreps;

MOMapping	moMappingNoPTRef;		// MOMapping (internal part without PTRefs)

private:
	ConfigurationSet	genRefSpace(ConfigurationSet, INT maxExc);
};



inline
Fort31File	RefSelInput::getMOIntegralFile() const
{	return MOIntegralFile;	}


inline
INT	RefSelInput::getFirstGuessConfs() const
{	return firstGuessConfs;	}

	

inline
double	RefSelInput::getReferenceThreshold() const
{	return ReferenceThreshold;	}



inline
const MRMOs	*RefSelInput::getMRMOs() const
{	return mrmos;	}




inline
const MOEquivalence	*RefSelInput::getMOEquivalence() const
{	return moEquivalence;	}


#endif
