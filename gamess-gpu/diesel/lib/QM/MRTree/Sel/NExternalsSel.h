//***********************************************************************
//
//	Name:			NExternalsSel.h
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


#ifndef __NExternalsSel_H
#define __NExternalsSel_H

#include "../../../../config.h"


#include "../Base/NExternalsBase.h"
#include "../Base/MRTreeBase.h"
#include "../Container/IndexedContainerTree.h"

#include "../Container/TreeRoot.h"

#include "../../Configuration/DiffConf.h"
#include "../../MO/MOType.h"
#include "../../MO/MRMOs.h"
#include "../../Math/SpinEigenFunctionDegeneration.h"

#include "../../Configuration/ConfigurationSet.h"


class	String;

template <class KeyType, class ObjectType> class	VarSizeReadOnlyCache;

template <class TMOType> class	Configuration;
template <class TMOType> class	TableCase;
class	TupelStructureSel;
class	InternalConfsSel;

class	BinomialCoefficient;

class	MOEquivalence;
class extEntry;


class NExternalsSel : 
	virtual public NExternalsBase<InternalConfsSel>,
	virtual public IndexedContainerTree<void, InternalConfsSel>,
	public MRTreeBase<void, InternalConfsSel> {
public:
	NExternalsSel(
		const MRMOs *mrmos,
		INT NumberOfElectrons,
		INT Multiplicity,
		INT ExcitationLevel,
		INT NumberOfRoots,
		const ConfigurationSet &confSet);

	NExternalsSel(istream &s, const MOIrReps *moIrReps);
	~NExternalsSel();
	
	//----------------------------------------------------------------	

	void	addExtEntry(const Configuration<MOType> &, const RootEnergies &);
	const extEntry *	getExtEntry(Configuration<MOType> conf) const;

	//----------------------------------------------------------------	

	void	cutTreeLevel2();

	void	threshCut(EnergyType Threshold, INT divideByCSFs);

	//----------------------------------------------------------------	

	void	writeToStream(ostream & s) const;

	friend ostream& operator<<(ostream & s, const NExternalsSel &);
	friend istream& operator>>(istream & s, NExternalsSel &);

	//----------------------------------------------------------------	

private:
};




#endif
