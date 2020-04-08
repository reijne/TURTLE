//***********************************************************************
//
//	Name:			MRConfInput.h
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

#ifndef __MRCONFINPUT_H
#define __MRCONFINPUT_H

#include <iostream>

#include "../../../config.h"

#include "../../../lib/Container/String.h"

#include "../../../lib/Container/SLList.h"
#include "../../../lib/QM/RepresentationMatrices/DataTypes.h"
#include "../../../lib/QM/MRTree/EnergyType.h"
#include "../../../lib/QM/Configuration/Configuration.h"
#include "../../../lib/QM/IO/Fortran/Fort31File.h"
#include "../../../lib/QM/Configuration/ConfigurationSet.h"
#include "EnergyMap.h"
#include "../../../lib/QM/MO/MOMapping.h"

class MRMOs;
class MORestriction;
class MOEquivalence;
class MOStatistics;

class MRConfInput {
public:
	MRConfInput();
	~MRConfInput();

friend class yyFlexLexer;


	// copy-constructor
	MRConfInput(const MRConfInput &);

	MRConfInput & operator = (const MRConfInput &);

	//--------------------------------------------------------------------------

	void	initInternalExternal();

	//--------------------------------------------------------------------------
	Fort31File	getMOIntegralFile() const;

	INT	getNumberOfElectrons() const;

	INT	getExcitationLevel() const;

	INT	getMultiplicity() const;
	
	INT	getFirstGuessConfs() const;
	
	INT	getSelectInternal() const;
	
	INT	getSelectNthExcitation() const;
	
	INT	getStorePTEnergy() const;
	INT	getStorePTCoef() const;

	IrRep	getIrRep() const;
	
	double getActiveReferenceThreshold() const;

	INT	getNumberOfSelectionThresholds() const;
	EnergyType	getSelectionThreshold(INT i) const;
	const char	*getSelectionThresholdString(INT i) const;
	const EnergyType	*getSelectionThresholdP() const;

	INT getRandThresh() const;
	double getRandomProb() const;

	ConfigurationSet &	getRefConfSet();
	ConfigurationSet &	getPTRefConfSet();

//	const Configuration<MOType>	*getRefConfs() const;
//	void	setRefConfs(const ConfigurationSet &confSet);

	const MRMOs	*getMRMOs() const;


	const MOStatistics	*getMOStatistics() const;

	INT	getNumberOfRoots() const;
	INT getRootNumber(INT i) const;

	RootType getRootEnergy(INT i) const;
	const RootType* getRootEnergyP() const;

	void	setVerbosity(const char *);

	const MORestriction	*getMORestriction() const;
	void	setMORestriction(const char *);

	const MOEquivalence	*getMOEquivalence() const;
	void	setMOEquivalence(const char *, INT auto = 0);

	MOMapping &	getMOMappingNoPTRef();

	EnergyMap::EstimationMode	getEstimationMode() const;
	//--------------------------------------------------------------------------

	friend ostream& operator<<(ostream & s, const MRConfInput &);
	friend istream& operator>>(istream & s, MRConfInput &);

	//--------------------------------------------------------------------------

protected:
Fort31File	MOIntegralFile;			// MO integral filename
INT	NumberOfElectrons;

INT	ExcitationLevel;				// number of processed excitations

INT	selectInternal;					// flag if always to select internal space

INT	selectNthExcitation;			// bit field coding which excitation levels
									// are selected unconditionally

INT	storePTEnergy;					// flag if PT energies are stored in
									// configuration tree

INT	storePTCoef;					// flag if PT coefs are stored in
									// configuration tree


MOStatistics	*moStatistics;		// energy statistics of MO pairs

INT	firstGuessConfs;				// number of configurations used in
									// first guess minus #roots

INT	autoRef;						// determine reference configurations
									// automatically
									
INT	autoEquiv;						// determine degenerated MOs
									// automatically
									
INT	Multiplicity;					// multiplicity of state

IrRep	irrep;						// number of irreducible representation

EnergyType	*SelectionThreshold;	// threshold for selection procedure
char	**SelectionThresholdString;	// same, but as a character string
INT	nSelectionThresholds;			// number of thresholds


double	randomProb;						// random selection threshold


ConfigurationSet	Refs;			// reference configurations
ConfigurationSet	PTRefs;			// PT reference configurations

INT	NumberOfRoots;
INT	NumberOfRootsE;
INT	*rootNumbers;
RootType	*rootEnergies;
MRMOs	*mrmos;


SLList<INT>	annihilatorSpace;
SLList<INT>	creatorSpace;
INT	activeSpaceExcitationLevel;
INT	maxRefOpenShells;
double	activeReferenceThreshold;

MORestriction	*moRestriction;
MOEquivalence	*moEquivalence;

EnergyMap::EstimationMode	estimationMode;

MOMapping	moMappingNoPTRef;		// MOMapping (internal part without PTRefs)

private:
	INT	checkConfs(ConfigurationSet &,
		INT &els, INT &elsErrors, INT &errors);
	void	genActiveSpace(
		Configuration<MOType> conf, INT *active,
		SLList<INT>, SLList<INT>);
	ConfigurationSet	genRefSpace(ConfigurationSet, INT maxExc);
	INT	isAutoRef() const;
	INT	isAutoEquiv() const;
};



inline
Fort31File	MRConfInput::getMOIntegralFile() const
{	return MOIntegralFile;	}

inline
INT	MRConfInput::getNumberOfElectrons() const
{	return NumberOfElectrons;	}

inline
INT	MRConfInput::getExcitationLevel() const
{	return ExcitationLevel;	}

inline
INT	MRConfInput::getMultiplicity() const
{	return Multiplicity;	}

inline
INT	MRConfInput::getFirstGuessConfs() const
{	return firstGuessConfs;	}

	
inline
INT	MRConfInput::isAutoRef() const
{	return autoRef;	}
	
inline
INT	MRConfInput::isAutoEquiv() const
{	return autoEquiv;	}

inline
INT	MRConfInput::getSelectInternal() const
{	return	selectInternal;	}
	
inline
INT	MRConfInput::getSelectNthExcitation() const
{	return	selectNthExcitation;	}
	
inline
INT	MRConfInput::getStorePTEnergy() const
{	return	storePTEnergy;	}

inline
INT	MRConfInput::getStorePTCoef() const
{	return	storePTCoef;	}

inline
double MRConfInput::getActiveReferenceThreshold() const
{	return activeReferenceThreshold;	}

inline
INT MRConfInput::getRandThresh() const
{	return static_cast<INT>((1-randomProb)*RAND_MAX);	}

inline
double MRConfInput::getRandomProb() const
{	return randomProb;	}

inline
IrRep	MRConfInput::getIrRep() const
{	return irrep;	}

inline
INT	MRConfInput::getNumberOfSelectionThresholds() const
{	return nSelectionThresholds;	}

inline
EnergyType	MRConfInput::getSelectionThreshold(INT i) const
{	return SelectionThreshold[i];	}

inline
const char	*MRConfInput::getSelectionThresholdString(INT i) const
{	return SelectionThresholdString[i];	}

inline
const EnergyType	*MRConfInput::getSelectionThresholdP() const
{	return SelectionThreshold;	}


inline
ConfigurationSet &	MRConfInput::getRefConfSet()
{	return Refs;	}

inline
ConfigurationSet &	MRConfInput::getPTRefConfSet()
{	return PTRefs;	}


inline
const MORestriction	*MRConfInput::getMORestriction() const
{	return moRestriction;	}

inline
const MOEquivalence	*MRConfInput::getMOEquivalence() const
{	return moEquivalence;	}

inline
const MOStatistics	*MRConfInput::getMOStatistics() const
{	return moStatistics;	}

inline
EnergyMap::EstimationMode	MRConfInput::getEstimationMode() const
{	return	estimationMode;	}

/*
inline
const Configuration<MOType>	*MRConfInput::getRefConfs() const
{	return RefConfs;	}
*/

inline
const MRMOs	*MRConfInput::getMRMOs() const
{	return mrmos;	}

inline
INT	MRConfInput::getNumberOfRoots() const
{	return NumberOfRoots;	}

inline
INT	MRConfInput::getRootNumber(INT i) const
{	return rootNumbers[i];	}

inline
RootType MRConfInput::getRootEnergy(INT i) const
{	return rootEnergies[i];	}

inline
const RootType* MRConfInput::getRootEnergyP() const
{	return rootEnergies;	}


inline
MOMapping &	MRConfInput::getMOMappingNoPTRef()
{	return	moMappingNoPTRef;	}


#endif
