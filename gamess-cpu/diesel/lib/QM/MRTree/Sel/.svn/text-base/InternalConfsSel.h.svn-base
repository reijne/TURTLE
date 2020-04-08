//***********************************************************************
//
//	Name:			InternalConfsSel.h
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


#ifndef __InternalConfsSel_H
#define __InternalConfsSel_H

#include "../../../../config.h"


#include "../Base/MRTreeBase.h"
#include "../Base/InternalConfsBase.h"

#include "TupelStructureSel.h"

#include "../Container/SetContainerTree.h"

template <class TMOType> class	Configuration;
class NExternalsSel;
class	MOEquivalence;

class InternalConfsSel : 
	public InternalConfsBase,
	public SetContainerTree<NExternalsSel, TupelStructureSel>,
	public MRTreeBase<NExternalsSel, TupelStructureSel> {
public:
	InternalConfsSel(INT nExt, NExternalsSel *parent) :
                Tree<NExternalsSel>(parent),
		InternalConfsBase(nExt) {   generated = 0;}
//	InternalConfsSel(NExternalsSel *parent);
	InternalConfsSel(istream &s);
	~InternalConfsSel();
	
	//----------------------------------------------------------------	

	void	addExtEntry(const Configuration<MOType> &, const RootEnergies &);
	const extEntry *	getExtEntry(Configuration<MOType> conf) const;
	 
	//----------------------------------------------------------------	


	INT	getGenerated() const;
	INT	addToGenerated(INT n);

	//----------------------------------------------------------------	

	void	cutTreeLevel2();

	void	threshCut(EnergyType Threshold, INT divideByCSFs);

	//----------------------------------------------------------------	

	void	writeToStream(ostream & s) const;

	friend ostream& operator<<(ostream & s, const InternalConfsSel &);
	friend istream& operator>>(istream & s, InternalConfsSel &);

	//----------------------------------------------------------------	

private:
INT	generated;		//	number of generated configurations
};


inline
INT	InternalConfsSel::getGenerated() const
{	return generated;	}

inline
INT	InternalConfsSel::addToGenerated(INT n)
{	return generated += n;	}


#endif
