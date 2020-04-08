//***********************************************************************
//
//	Name:			Selector.h
//
//	Description:	
//
//	Author:			Michael Hanrath
//
//	Version:		0.0
//
//	Date:			23.03.1997
//
//
//
//
//
//***********************************************************************



#ifndef __Selector_h
#define __Selector_h

#include "../../../config.h"

#include "../../../lib/QM/Configuration/DiffConf.h"
#include "../../../lib/QM/Configuration/TableCase.h"

#include "../../../lib/QM/MRCIMatrix/CICalculation.h"

#include "../../../lib/QM/RepresentationMatrices/CIVectors.h"


#include "MRConfInput.h"

#include "../../../lib/QM/MRTree/Sel/NExternalsSel.h"
#include "../../../lib/QM/MRTree/Set/NExternalsSet.h"
#include "../../../lib/QM/MO/Iterators/MOList.h"


template <class Key, class ContainedObject> class VarSizeReadOnlyCache;


template <class MatrixType> class RepresentationMatrixFortranInterface;
class TupelStructureSel;
template <class MatrixType, class VectorType> class PTReferenceSpace;
template <class MatrixType, class VectorType> class EnlargeReferenceSpace;

template <class T> class Histogram;

template <class MatrixType, class VectorType>
class Selector :
	public CICalculation<MatrixType, VectorType>
 {
public:
	Selector(MRConfInput mrinp, NExternalsSet *preSelConfs = NULL,
		INT cacheEntries = 2048, INT cacheMemory = (1 << 20));
	~Selector();



	void	useAnnihilators();
	void	useCreators(
		const PTReferenceSpace<MatrixType, VectorType> *PTReference,
		EnergyType Threshold);


	NExternalsSel	*getGenConfs();

	EnergyType	getPreESum(INT root) const {	return preESum[root];	}

	EnergyType	calcSingle(Configuration<MOType>, Configuration<MOType>);

	//----------------------------------------------------------------	

	void printResults() const;

	VarSizeReadOnlyCache<TableKey, HMatElements<MatrixType, VectorType> >	*getCache() const;

private:
	INT	checkAmbiguity(INT n, INT *isInternal, IrRep *irrep) const;
	void select(
		TupelStructureSel *tupel,
		const MOList & molist,
		EnergyType Threshold);
	

VarSizeReadOnlyCache<TableKey, HMatElements<MatrixType, VectorType> >	*cache;


NExternalsSel	*genConfs;					// pointer to tree of generated/
											// selected confs

const PTReferenceSpace<MatrixType, VectorType> *PTReference;		// pointer to pertubation theory
											// reference space


EnlargeReferenceSpace<MatrixType, VectorType> *enlargeRef;			// energy estimation by reference space
											// enlargement


MRConfInput	mrinp;							// user input data

TableCase<GeneralizedMO>	*tablecaseGen;	// table case based on generalized MOs

DiffConf<MOType>	dc;						// difference configuration
DiffConf<GeneralizedMO>	dcGen;				// difference configuration based on 
											// generalized MOs


RepresentationMatrixFortranInterface<MatrixType>	*repMatFInt;

HMatElements<MatrixType, VectorType>	**repMatsP5;		// pointer to array of P=5
											// representation matrices for
											// a different number of open shells
											// (used in calculation of 
											// diagonal elements)

INT	minOpenShells;							// minimum number of open shells
											// for given multiplicity
											
INT	*genConfCount;							// array of numbers of generated
											// configurations with n externals

INT	*genRestrictedConfCount;				// dito. restricted

INT	*genCSFCount;							// array of numbers of generated
											// configurations with n externals											

INT	*genRestrictedCSFCount;					// dito. restricted

EnergyType	*RestrictedESum;				// energy sum from restricted configurations



NExternalsSet	*preSelConfs;				// if NULL: perform selection
											// if not NULL: recalculation mode
											// (pointer to tree of 
											// preselected confs)
											
EnergyType	*preESum;						// energy sum in recalculation mode

Histogram<EnergyType>	*diagHist;			// histogram of diagonal elements
};


template <class MatrixType, class VectorType>
inline
NExternalsSel *Selector<MatrixType, VectorType>::getGenConfs()
{	return genConfs;	}

template <class MatrixType, class VectorType>
inline
VarSizeReadOnlyCache<TableKey, HMatElements<MatrixType, VectorType> >	*Selector<MatrixType, VectorType>::getCache() const
{	return cache;	}


#endif
