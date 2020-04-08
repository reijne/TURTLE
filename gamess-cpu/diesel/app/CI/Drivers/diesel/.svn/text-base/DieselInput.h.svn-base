//***********************************************************************
//
//	Name:			dieselInput.h
//
//	Description:	stores input for diesel driver 
//
//	Author:			Michael Hanrath
//
//	Version:		0.0
//
//	Date:			09.11.1998
//
//
//
//
//
//***********************************************************************

#ifndef __dieselInput_H
#define __dieselInput_H

#include "../../../../config.h"

#include <iostream>
#include <vector>
using std::vector;

#include "../../../../lib/Container/String.h"

#include "../../../../lib/QM/RepresentationMatrices/DataTypes.h"

#include "../../../../lib/QM/MRTree/EnergyType.h"
#include "../../../../lib/QM/IO/Fortran/Fort31File.h"
#include "../../Selector/EnergyMap.h"

#include "../../../../lib/QM/MO/MOMapping.h"

#include "../../../../lib/QM/MRCIMatrix/MRMPH0Matrix.h"
#include "../../../../lib/QM/Configuration/ConfigurationSet.h"
#include "../../../../lib/QM/Symmetry/DegenPointgroups/RefSelMode.h"

#include "../../../../lib/Container/SLList.h"


class MRMOs;
class MORestriction;
class MOEquivalence;
class ConfigurationSet;


class DieselInput {
public:
	DieselInput();
	~DieselInput();

friend class yyFlexLexer;

typedef double	VectorType;


	//--------------------------------------------------------------------------

	Fort31File	getMOIntegralFile() const;

	ConfigurationSet &	getRefConfSet();


	const MRMOs	*getMRMOs() const;

	const SLList<String> &	getSelectionThresholds() const;
	const vector<INT> &	getRoots(INT irrep) const;
	const SLList<INT> &	getIrreps() const;
	const SLList<INT> &	getMultiplicities() const;
	const SLList<INT> &	getNthExcitation() const;


	const MORestriction	*getMORestriction() const;
	void	setMORestriction(const char *);

	const MOEquivalence	*getMOEquivalence() const;
	void	setMOEquivalence(const char *, INT auto = 0);


	void	setVerbosity(const char *);


	void	start();
	
	
	//--------------------------------------------------------------------------

	friend ostream& operator<<(ostream & s, const DieselInput &);
	friend istream& operator>>(istream & s, DieselInput &);

	//--------------------------------------------------------------------------

protected:
String	MOLCASRootDir;				// MOLCAS base directory
String	TempDir;					// temporary directory
Fort31File	MOIntegralFile;			// MO integral filename
INT	NumberOfElectrons;

INT	ExcitationLevel;				// number of processed excitations

INT	selectInternal;					// flag if always to select internal space

INT	autoRef;						// determine reference configurations
									// automatically


RefSelMode::Mode	refSelMode;
								
INT	autoEquiv;						// determine degenerated MOs
									// automatically

INT	nProc;							// number of processors to be used

INT MaxIters;								// maximum number of iterations
//LONG_INT	MaxHamiltonStorageMem;					// maximum memory consumption 
INT      MaxHamiltonStorageMem;                                  // maximum memory consumption


INT	maxRefGenIters;

INT	firstGuessConfs;				// number of configurations used in
									// first guess minus #roots

EnergyType	ConvergenceEnergyChange;		// convergence threshold (energy)

VectorType	ConvergenceEigenvectorChange;	// convergence threshold (eigen vector)

MORestriction	*moRestriction;
MOEquivalence	*moEquivalence;

EnergyMap::EstimationMode	estimationMode;

vector<INT>	rootsDefault;
vector<vector<INT> >	roots;
INT	preScanRoots;
SLList<INT>	irreps;
SLList<INT>	multiplicities;
SLList<String>	thresholds;
SLList<INT>	nthExcitation;

SLList<String>	fullMRCIExtrapolation;
SLList<String>	propertyThresholds;
String	orbitalFile;
SLList<String>	store;


ConfigurationSet	Refs;			// reference configurations

INT	rootHoming;						// flag if root homing is active


INT	usePreviousSelection;			// flag if already performed selection
									// should be used 
									// (useful with potential surfaces, gradients)

INT homogeneousSelection;			// flag if "selsym" is applied to
									// selected configurations

SLList<INT>	annihilatorSpace;
SLList<INT>	creatorSpace;
INT	activeSpaceExcitationLevel;
INT	maxRefOpenShells;
double	activeReferenceThreshold;

INT	useNaturalOrbitals;
INT	averagedNaturalOrbitals;
String	NaturalOrbitalSelectionThreshold;

String	MRPTInhomogenityThreshold;
SLList<String>	MRPTthresholds;

VectorType	RefThreshold;				// threshold for reference selection
VectorType	PTRefThreshold;				// threshold for PT reference selection
VectorType	NatOrbRefThreshold;			// threshold for natural orbital 
										// generation reference selection

MRMOs	*mrmos;

MRMPH0Matrix<double, double>::ProjectionMode projectionMode;

private:
	// copy-constructor
	DieselInput(const DieselInput &);

	DieselInput & operator = (const DieselInput &);

	void	createReferenceSpace(INT multiplicity);
	void	createNatOrbs(Pix iMult, Pix iIrrep, INT &nNatOrbCalcs);

	void	doDIESEL(
		INT multiplicity,
		INT irrep,
		const SLList<String> &thresholds,
		INT doNatOrb);

	void	writeRefGuessInput(
		String filename,
		vector<INT>	roots,
		INT multiplicity) const;

	void	writeRefSelInput(
		String filename, SLList<INT> irreps) const;

	void	writeSelectorInput(
		String filename,
		vector<INT>	roots,
		SLList<INT>	nthExcitation,
		INT multiplicity,
		INT irrep,
		const SLList<String> &thresholds,
		const ConfigurationSet &,
		const ConfigurationSet &,
		bool supressActive = false) const;

	void	writeMRPTInput(
		String filename,
		vector<INT>	roots,
		const SLList<String> &thresholds,
		INT mp3) const;
		
	void	writeDiagonalisatorInput(
		String filename,
		VectorType	RefThreshold,
		VectorType	PTRefThreshold,
		vector<INT>	roots
		) const;
};


inline
ConfigurationSet &	DieselInput::getRefConfSet()
{	return Refs;	}

inline
Fort31File	DieselInput::getMOIntegralFile() const
{	return MOIntegralFile;	}

inline
const SLList<String> &	DieselInput::getSelectionThresholds() const
{	return thresholds;	}

inline
const vector<INT> &	DieselInput::getRoots(INT irrep) const
{	return roots[irrep];	}

inline
const SLList<INT> &	DieselInput::getIrreps() const
{	return irreps;	}

inline
const SLList<INT> &	DieselInput::getMultiplicities() const
{	return multiplicities;	}

inline
const SLList<INT> &	DieselInput::getNthExcitation() const
{	return nthExcitation;	}


inline
const MORestriction	*DieselInput::getMORestriction() const
{	return moRestriction;	}

inline
const MOEquivalence	*DieselInput::getMOEquivalence() const
{	return moEquivalence;	}



inline
const MRMOs	*DieselInput::getMRMOs() const
{	return mrmos;	}


#endif
