//***********************************************************************
//
//	Name:			MOIrReps.h
//
//	Description:	
//
//	Author:			Michael Hanrath
//
//	Version:		0.0
//
//	Date:			12.08.1996
//
//
//
//
//
//***********************************************************************

#ifndef __MOIRREPS_H
#define __MOIRREPS_H

#include "../../../config.h"

#include "../Symmetry/IrRep.h"
#include "MOType.h"
#include "../IO/Fortran/Fort31File.h"
#include "RISym.h"

#include <iostream>
using std::istream;
using std::ostream;

class MOIrReps {
public:
	MOIrReps();
	MOIrReps(istream &);
	MOIrReps(Fort31File f31);
	~MOIrReps();


	// copy-constructor
	MOIrReps(const MOIrReps &);

	// Zuweisungs-Operator
	MOIrReps & operator = (const MOIrReps &);

	
	//----------------------------------------------------------------	


	MOType	getMaxMO() const;
	
	INT	getNumberOfIrReps() const;
	
	INT	getInIrRep(IrRep i) const;

	MOType	getStartMO(IrRep i) const;
	MOType	getEndMO(IrRep i) const;
	
	IrRep	getIrRep(MOType mo) const;
	IrRep	getIrRepP1(MOType mo) const;
	
	
	MOType	getNumberInIrRep(MOType mo) const;


	//----------------------------------------------------------------	

	const INT *	getInIrRepP() const;
	const IrRep *	getMOSymmetryP() const;
	const IrRep *	getProdTabP() const;
	const MOType *	getNumberInIrRepP() const;
	const MOType *	getStartMOP() const;
	
	//----------------------------------------------------------------	

	IrRep	getProd(IrRep, IrRep) const;
//	IrRep &	getProd(IrRep, IrRep);

	//----------------------------------------------------------------

	friend ostream& operator<<(ostream & s, const MOIrReps &);

	ostream & writeToStream(ostream &);

	//----------------------------------------------------------------
	
protected:
	void	init();

static const INT maxIrreps;
static const IrRep prodTab[];


MOType	maxMO;			// maximum MO number
INT	IrReps;				// number of IrReps
INT	*inIrRep;			// number of MOs in irrep
IrRep	*MOSymmetry;	// translation table MO -> irrep
MOType	*MOinIrRep;		// translation table MO -> MO in irrep
MOType	*IrRepStart;	// start MO of Irrep
//IrRep	*prodTab;		// product table
};




inline
IrRep	MOIrReps::getIrRep(MOType mo) const
{	return	MOSymmetry[mo];	}

inline
IrRep	MOIrReps::getIrRepP1(MOType mo) const
{	return	MOSymmetry[mo] + 1;	}

inline
INT	MOIrReps::getNumberOfIrReps() const
{	return	IrReps;	}

inline
INT	MOIrReps::getInIrRep(IrRep i) const
{	return inIrRep[i];	}

inline
MOType	MOIrReps::getMaxMO() const
{	return	maxMO;	}

inline
IrRep	MOIrReps::getProd(IrRep irrep1, IrRep irrep2) const
{	return	prodTab[irrep1*maxIrreps + irrep2];	}

/*inline
IrRep &	MOIrReps::getProd(IrRep irrep1, IrRep irrep2)
{	return	prodTab[irrep1*maxIrreps + irrep2];	}
*/
inline
const INT *	MOIrReps::getInIrRepP() const
{	return inIrRep;	}

inline
const IrRep *	MOIrReps::getMOSymmetryP() const
{	return MOSymmetry;	}

inline
const MOType *	MOIrReps::getNumberInIrRepP() const
{	return	MOinIrRep;	}

inline
const IrRep *	MOIrReps::getProdTabP() const
{	return prodTab;	}

inline
const MOType *	MOIrReps::getStartMOP() const
{	return IrRepStart;	}

inline
MOType	MOIrReps::getNumberInIrRep(MOType mo) const
{	return MOinIrRep[mo-1];	}

inline
MOType	MOIrReps::getStartMO(IrRep i) const
{	return IrRepStart[i];	}

inline
MOType	MOIrReps::getEndMO(IrRep i) const
{	return IrRepStart[i]+inIrRep[i]-1;	}

#endif
