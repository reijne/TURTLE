//***********************************************************************
//
//	Name:			NExternalsSet.h
//
//	Description:	
//
//	Author:			Michael Hanrath
//
//	Version:		0.0
//
//	Date:			30.05.1997
//
//
//
//
//
//***********************************************************************


#ifndef __NExternalsSet_H
#define __NExternalsSet_H

#include "../../../../config.h"


#include "../Base/NExternalsBase.h"
#include "../Base/MRTreeBase.h"
#include "../Container/SetContainerTree.h"

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
class	TupelStructureSet;
class	InternalConfsSet;

class	BinomialCoefficient;

class NExternalsSet : 
	virtual public NExternalsBase<InternalConfsSet>,
	virtual public SetContainerTree<void, InternalConfsSet>,
	public MRTreeBase<void, InternalConfsSet> {
public:
	NExternalsSet();
	NExternalsSet(istream &s, const MRMOs *mrmos = NULL);
	NExternalsSet(
		const MRMOs *mrmos,
		INT NumberOfElectrons,
		INT Multiplicity,
		INT NumberOfRoots,
		const ConfigurationSet &refConfs);
	NExternalsSet(
		const MRMOs *mrmos,
		INT NumberOfElectrons,
		INT Multiplicity,
		INT NumberOfRoots,
		IrRep totalSymmetry);

	~NExternalsSet();
	

	// "virtual" constructor:
//	virtual NExternalsSet * new_Set();

	void	copy(NExternalsSet *, INT deep = 0);
	virtual SetContainerTree<void, InternalConfsSet> *	clone(INT deep = 0);
	
	//----------------------------------------------------------------	

	void	add(const Configuration<MOType> &);

	void	add(InternalConfsSet *);

	//----------------------------------------------------------------	

	INT checkSameRefs(const NExternalsSet & ext) const;
	
	void	cutTreeLevel2();

	//----------------------------------------------------------------
	
	void	writeToStream(ostream & s) const;

	friend ostream& operator<<(ostream & s, const NExternalsSet &);
	friend istream& operator>>(istream & s, NExternalsSet &);

	//----------------------------------------------------------------	

private:
};


/*
inline 
NExternalsSet *	NExternalsSet::
	new_Set()
{	return new NExternalsSet();	}
*/

inline 
void	NExternalsSet::copy(NExternalsSet *c, INT deep)
{
	*this = *c;
	if ( deep )
		SetContainerTree<void, InternalConfsSet>::clone(1);
}

inline
void	NExternalsSet::add(InternalConfsSet *t)
{
	SetContainerTree<void, InternalConfsSet>::add(t);
}

inline 
SetContainerTree<void, InternalConfsSet> *	NExternalsSet::clone(INT deep)
{	
	NExternalsSet * c = new NExternalsSet();
	c->copy(this, deep);
	return c;
}




#endif
