//***********************************************************************
//
//	Name:			NExternalsBase.h
//
//	Description:	
//
//	Author:			Michael Hanrath
//
//	Version:		0.0
//
//	Date:			09.03.1997
//
//
//
//
//
//***********************************************************************


#ifndef __NExternalsBase_h
#define __NExternalsBase_h

#include "../../../../config.h"


#include "../../Configuration/DiffConf.h"
#include "../../MO/MOType.h"
#include "../../MO/MRMOs.h"
#include "../../Math/SpinEigenFunctionDegeneration.h"

#include "../Container/Container.h"
#include "InternalConfsBase.h"
#include "../../Configuration/ConfigurationSet.h"

using std::ostream;
using std::istream;

class	String;

/*
template <class KeyType, class ObjectType> class	VarSizeReadOnlyCache;

template <class TMOType> class	Configuration;
template <class TMOType> class	TableCase;
class	TupelStructure;
class	InternalConfs;

class	BinomialCoefficient;
*/

template <class ContainedObjectType> class NExternalsBase;
template <class ContainedObjectType> 
         ostream& operator << (ostream & s, const NExternalsBase<ContainedObjectType> &);
template <class ContainedObjectType> 
         istream& operator >> (istream & s, NExternalsBase<ContainedObjectType> &);


template <class ContainedObjectType>
class NExternalsBase : virtual public Container<ContainedObjectType> {
public:
	NExternalsBase();
	NExternalsBase(
		const MOIrReps *moIrReps,
		INT NumberOfElectrons, INT Multiplicity);
	NExternalsBase(
		const MRMOs *mrmos,
		INT NumberOfElectrons, INT Multiplicity);
	NExternalsBase(istream &s, const MOIrReps *moIrReps);
	virtual ~NExternalsBase();
	
	// copy-constructor
	NExternalsBase(const NExternalsBase<ContainedObjectType> &);
	NExternalsBase & operator = (const NExternalsBase<ContainedObjectType> &);

	
	//----------------------------------------------------------------	
	
	INT	getMultiplicity() const;
	INT	getNumberOfElectrons() const;
	IrRep	getTotalSymmetry() const;

	//	get spin eigenfunction degeneration
	INT	getNumberOfSpinAdaptedFunctions(INT) const;

	//----------------------------------------------------------------	

	INT	getNumberOfTotalSpinAdaptedFunctions() const;
	INT	getNumberOfRefConfSpinAdaptedFunctions() const;
	INT	getNumberOfReferences() const;

	INT getNumberOfRoots() const;
	void setNumberOfRoots(INT n);
	INT	getRootNumber(INT i) const;
	void setRootNumber(INT i, INT n) const;
	const INT * getRootNumbersP() const;

	void	addRefSAFs(INT i);
	
	//----------------------------------------------------------------	

	const MRMOs *	getMRMOs() const;

	//----------------------------------------------------------------	

	void	calcInternalIntersection();
	const Configuration<MOType> & getInternalIntersection() const;
	ConfigurationSet	getRefConfSet() const;

	//----------------------------------------------------------------	

	void	writeToStream(ostream & s) const;
	
	friend ostream& operator << <ContainedObjectType> (ostream & s, const NExternalsBase<ContainedObjectType> &);
	friend istream& operator >> <ContainedObjectType> (istream & s, NExternalsBase<ContainedObjectType> &);

	//----------------------------------------------------------------	

protected:
	void	initAll();


INT	SAFs;				// total number spin adapted functions
INT	RefSAFs;			// number spin adapted functions from reference
						// configurations

INT	NumberOfRoots;		// number of roots the selection is based on
INT	*rootNumbers;		// array with number of i-th root


IrRep	totalSymmetry;
INT	electrons;			// number of electrons
INT	multiplicity;		// multiplicity
INT	References;			// number of references

MRMOs	*mrmos;			// information about MOs
SpinEigenFunctionDegeneration	*eigFuncs;
static const INT	maxSym = 8;


Configuration<MOType>	internalIntersection;
						// internal occupation pattern common to all
						// configurations in tree
};

template <class ContainedObjectType>	ostream& operator << (ostream & s, const NExternalsBase<ContainedObjectType> &);
template <class ContainedObjectType>	istream& operator >> (istream & s, NExternalsBase<ContainedObjectType> &);

template <class ContainedObjectType>
inline
const MRMOs *	NExternalsBase<ContainedObjectType>::getMRMOs() const
{	return	mrmos;	}

template <class ContainedObjectType>
inline
void	NExternalsBase<ContainedObjectType>::addRefSAFs(INT i)
{	RefSAFs += i;	}

template <class ContainedObjectType>
inline
INT	NExternalsBase<ContainedObjectType>::getNumberOfTotalSpinAdaptedFunctions() const
{	return SAFs;	}

template <class ContainedObjectType>
inline
INT	NExternalsBase<ContainedObjectType>::getNumberOfRefConfSpinAdaptedFunctions() const
{	return RefSAFs;	}

template <class ContainedObjectType>
inline
INT	NExternalsBase<ContainedObjectType>::getNumberOfSpinAdaptedFunctions(
	INT openShells) const
{	return (*eigFuncs)(openShells);	}

template <class ContainedObjectType>
inline
INT	NExternalsBase<ContainedObjectType>::getNumberOfReferences() const
{	return References;	}

template <class ContainedObjectType>
inline
INT	NExternalsBase<ContainedObjectType>::getMultiplicity() const
{	return	multiplicity;	}

template <class ContainedObjectType>
inline
INT	NExternalsBase<ContainedObjectType>::getNumberOfElectrons() const
{	return	electrons;	}

template <class ContainedObjectType>
inline
IrRep	NExternalsBase<ContainedObjectType>::getTotalSymmetry() const
{	return	totalSymmetry;	}


template <class ContainedObjectType>
inline
INT NExternalsBase<ContainedObjectType>::getNumberOfRoots() const
{	return	NumberOfRoots;	}

template <class ContainedObjectType>
inline
void NExternalsBase<ContainedObjectType>::setNumberOfRoots(INT n)
{
	if ( rootNumbers )
		delete rootNumbers;
	NumberOfRoots = n;
	rootNumbers = new INT[n];
}

template <class ContainedObjectType>
inline
INT	NExternalsBase<ContainedObjectType>::getRootNumber(INT i) const
{	return rootNumbers[i];	}

template <class ContainedObjectType>
inline
void	NExternalsBase<ContainedObjectType>::setRootNumber(INT i, INT n) const
{	rootNumbers[i] = n;	}

template <class ContainedObjectType>
inline
const INT * NExternalsBase<ContainedObjectType>::getRootNumbersP() const
{	return rootNumbers;	}

template <class ContainedObjectType>
inline
const Configuration<MOType> & NExternalsBase<ContainedObjectType>::getInternalIntersection() const
{	return internalIntersection;	}

#endif
