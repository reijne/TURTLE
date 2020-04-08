//***********************************************************************
//
//	Name:			NExternalsDiag.h
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


#ifndef __NExternalsDiag_H
#define __NExternalsDiag_H

#include "../../../../config.h"

#include "../Base/NExternalsBase.h"

#include "DiagTreeBase.h"
#include "../Container/TreeRoot.h"

#include "../../Configuration/DiffConf.h"
#include "../../MO/MOType.h"
#include "../../MO/MRMOs.h"
#include "../../Math/SpinEigenFunctionDegeneration.h"

//#include "MRConfInput.h"


class	String;
//class	CICalculation;

template <class KeyType, class ObjectType> class	VarSizeReadOnlyCache;

template <class TMOType> class	Configuration;
template <class TMOType> class	TableCase;
class	TupelStructureDiag;
class	InternalConfsDiag;

class	BinomialCoefficient;

class	NExternalsSet;

class NExternalsDiag : 
	virtual public NExternalsBase<InternalConfsDiag>,
	public DiagTreeBase<void, InternalConfsDiag> {
public:
	NExternalsDiag(
		const MOIrReps *moIrReps,
		INT NumberOfElectrons, INT Multiplicity,
		INT HighestOrder = 2) :
		NExternalsBase<InternalConfsDiag>(moIrReps, NumberOfElectrons, Multiplicity), 
                Tree<void>((void *) NULL),
		IndexedContainer<InternalConfsDiag>(HighestOrder+1) {}

	NExternalsDiag(istream &s, const MOIrReps *moIrReps);
/*	 :
		NExternalsBase(s), 
		IndexedContainer<InternalConfsDiag>(s),
		Tree<void>((void *) NULL) {}
*/	
//	NExternalsDiag(MRConfInput &mrconf);

	NExternalsDiag(const NExternalsDiag &);

	NExternalsDiag(const NExternalsSet &, INT referenceFlag = 0);

/*	NExternalsDiag(
		const ConfigurationSet &references, 
		const MRMOs *,
		INT nElectrons,
		INT Multiplicity,
		INT numberOfRoots
		);
*/

	~NExternalsDiag() {}

	
	//----------------------------------------------------------------	

	NExternalsDiag projectOnReferences() const;

	//----------------------------------------------------------------	



	//	Problem with GCC 2.8.1
	
/*	template <class MatrixAction, class MatrixElementCalculator> static void	twoDimIteration(
		const NExternalsDiag *, 
		const NExternalsDiag *,
		MatrixAction &, 
		const MatrixElementCalculator &dummy,
		INT highestExcitation,
		INT	upperTriangular);
*/
/*
	template <class MatrixAction, class MatrixElementCalculator> void	twoDimIteration(
		MatrixAction &action, 
		const MatrixElementCalculator &dummy,
		INT highestExcitation,
		INT	upperTriangular, 
		const NExternalsDiag *b) const
	{	twoDimIteration(this, b, action, dummy, highestExcitation, upperTriangular);	}

	template <class MatrixAction, class MatrixElementCalculator> void	twoDimIteration(
		MatrixAction &action, 
		const MatrixElementCalculator &dummy,
		INT highestExcitation,
		INT	upperTriangular) const
	{	twoDimIteration(this, this, action, dummy, highestExcitation, upperTriangular);	}
*/

	//----------------------------------------------------------------	
	
	void	writeToStream(ostream & s) const;

	friend ostream& operator<<(ostream & s, const NExternalsDiag &);
	friend istream& operator>>(istream & s, NExternalsDiag &);

	//----------------------------------------------------------------	

private:
	void	initAll();
	
};



#endif
