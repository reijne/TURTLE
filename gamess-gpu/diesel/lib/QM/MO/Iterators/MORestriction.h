//***********************************************************************
//
//	Name:			MORestriction.h
//
//	Description:	restrictions for MO occupation patterns
//
//	Author:			Michael Hanrath
//
//	Version:		0.0
//
//	Date:			08.03.1998
//
//
//
//
//
//***********************************************************************


#ifndef __MORestriction_h
#define __MORestriction_h

#include "../../../../config.h"

#include "../../../Container/BitSet.h"
#include "../../../Container/String.h"


#include "../MOType.h"

using std::ostream;

template <class TMOType> class	Configuration;
class MOIrReps;

class MORestriction {
friend class MORestrictionState;
public:
	MORestriction(INT maxMO, INT nRestricts);
	MORestriction(const Configuration<MOType> &conf, 
		const MOIrReps &moirreps,
		INT activeVirtual = 4);
	MORestriction(
		String restString);
	~MORestriction();

	MORestriction(const MORestriction &);

	void	setMaxMO(INT maxMO);

	INT	check(const Configuration<MOType> &conf) const;
	
	INT	isIllegal() const;
	
	String	getRestriction() const;

//	MORestriction &	operator |= (const MORestriction &);

	friend ostream & operator << (ostream &, const MORestriction &);


private:
	MORestriction & operator = (const MORestriction &);

	void	getOperatorOccupation(
		const char *p,
		INT	&op, INT &occupation) const;

	INT	getOperator(const char *) const;

	String	substRange(String, char separator = ',') const;
	

	
struct	TRestriction {
	BitSet	mask;					//	MO bit mask corresponding to 
									//	occupation array
	INT		occupation;				//	occupation number
	INT		op;						//	operator (<, <=, =, >, >=, !=)
	} *restrictions;
	
INT		maxMO;						//	highest MO number
INT		nRestrictions;				//	number of separate restrictions

static char *Operators[];
static const  INT	nOperators = 6;
};




#endif
