//***********************************************************************
//
//	Name:			TupelStructureDiag.h
//
//	Description:	
//
//	Author:			Michael Hanrath
//
//	Version:		0.1
//
//	Date:			21.06.1996
//					09.03.1997
//
//
//
//
//
//***********************************************************************


#ifndef __TupelStructureDiag_H
#define __TupelStructureDiag_H

#include "../../../../config.h"



#include "../Base/TupelStructureBase.h"
#include "../../Configuration/Configuration.h"
#include "extMOsDiag.h"
#include "DiagTreeBase.h"


class	extMOsDiag;
class	InternalConfsDiag;
class	TupelStructureSet;

class TupelStructureDiag : 
	public DiagTreeBase<InternalConfsDiag, extMOsDiag>,
	public TupelStructureBase {
//	public ContainerTree<TupelStructureDiag, extMOsDiag> {
public:
	TupelStructureDiag(INT n, INT symmetry, 
		INT NumberOfOpenShells, INT NumberOfClosedShells, MOType *pShells,
		InternalConfsDiag *parent, INT ReferenceFlag = 0) :
                Tree<InternalConfsDiag>(parent),
                IndexedContainer<extMOsDiag>(n),
                TupelStructureBase(
                symmetry, NumberOfOpenShells, NumberOfClosedShells,
                        pShells, ReferenceFlag) {}

	TupelStructureDiag(INT n, INT symmetry, 
		Configuration<MOType> conf,
		InternalConfsDiag *parent, INT ReferenceFlag = 0) :
                Tree<InternalConfsDiag>(parent),
                IndexedContainer<extMOsDiag>(n),
		TupelStructureBase(
		symmetry, conf, ReferenceFlag) {}

	TupelStructureDiag(const TupelStructureDiag &);
	
	TupelStructureDiag(istream &s);


	TupelStructureDiag(const TupelStructureSet &, InternalConfsDiag *, INT referenceFlag);

	~TupelStructureDiag() {}

	//----------------------------------------------------------------	

	INT getSAFStart() const;
	INT getSAFInc() const;
	
	INT	init(INT, const IrRep *);
		
	//----------------------------------------------------------------	

	void	writeToStream(ostream & s) const;

	friend ostream& operator<<(ostream & s, const TupelStructureDiag &);
	friend istream& operator>>(istream & s, TupelStructureDiag &);

	//----------------------------------------------------------------	

private:
INT	SAFStart;			// start of saf in CI-vector
INT	SAFInc;				// number of SAFs
};


inline
INT TupelStructureDiag::getSAFStart() const
{	return SAFStart;	}

inline
INT TupelStructureDiag::getSAFInc() const
{	return SAFInc;	}



#endif
