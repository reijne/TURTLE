//***********************************************************************
//
//	Name:			MOMapping.h
//
//	Description:	
//
//	Author:			Michael Hanrath
//
//	Version:		0.0
//
//	Date:			12.05.1997
//
//
//
//
//
//***********************************************************************

#ifndef __MOMapping_h
#define __MOMapping_h

#include "../../../config.h"

#include "MOType.h"
#include "GeneralizedMO.h"

#include <iostream>

class	NExternalsDiag;
template <class Type> class	CIVectors;



class MOMapping {
public:
	MOMapping();
	MOMapping(istream &s);
	MOMapping(INT maxMO, INT *internal, const MOIrReps &moirreps);
	~MOMapping();
	
	MOMapping(const MOMapping &);
	MOMapping & operator = (const MOMapping &);
	
	//----------------------------------------------------------------	
	
	INT	getMaxMO() const;
	
	MOType	getReal(MOType continuous) const;
	GeneralizedMO	getReal(GeneralizedMO continuous) const;
	MOType	getContinuous(MOType real) const;
	GeneralizedMO	getContinuous(GeneralizedMO real) const;
	
	//----------------------------------------------------------------	

	void	setInternal(MRMOs *) const;
	template <class Type>
	void	reorderNoHoles(NExternalsDiag * &, CIVectors<Type> * &,
		const MRMOs *) const;

	//----------------------------------------------------------------	

	void	writeToStream(ostream &s);

	friend ostream& operator<<(ostream & s, const MOMapping &);

	//----------------------------------------------------------------

private:
INT	maxMO;						// size of mapping tables
MOType	*mapRealToContinuous;		// mapping table Real --> Continuous
MOType	*mapContinuousToReal;		// mapping table Continuous --> Real
};


extern MOMapping moMapping;

inline
INT	MOMapping::getMaxMO() const
{	return maxMO;	}

inline
MOType	MOMapping::getReal(MOType continuous) const
{	
	if ( maxMO )
		return	mapContinuousToReal[continuous-1];
	else
		return continuous;
}

inline
GeneralizedMO	MOMapping::getReal(GeneralizedMO continuous) const
{	return continuous;	}

inline
MOType	MOMapping::getContinuous(MOType real) const
{
	if ( maxMO )
		return	mapRealToContinuous[real-1];
	else
		return real;
}

inline
GeneralizedMO	MOMapping::getContinuous(GeneralizedMO real) const
{	return real;	}





#endif
