//***********************************************************************
//
//	Name:			MOEquivalence.h
//
//	Description:	set of MOs to be handled equivalently
//
//	Author:			Michael Hanrath
//
//	Version:		0.0
//
//	Date:			10.04.1998
//
//
//
//
//
//***********************************************************************


#ifndef __MOEquivalence_h
#define __MOEquivalence_h

#include "../../../../config.h"

#include "../../../Container/String.h"
#include "../MOType.h"
#include "../../MRTree/EnergyType.h"

using std::ostream;

class ConfigurationSet;
template <class TMOType> class	Configuration;
class	MRMOs;

class MOEquivalence {
friend class MOEquivalenceProjected;
public:
	MOEquivalence(String equivString);
	MOEquivalence(const EnergyType *e, INT maxMO);
	~MOEquivalence();

	MOEquivalence(const MOEquivalence &);

struct TEquiv {
	INT n;							// number of equivalent MOs
	MOType	*mo;					// list of equivalent MOs
	};



	INT	getNumberOfEquivalences() const
	{	return	n;	}


	TEquiv getEquiv(INT i) const
	{	return equiv[i];	}

	String	getEquivalence() const;

	void	symmetrize(ConfigurationSet &, const MRMOs *);

	friend ostream & operator << (ostream &, const MOEquivalence &);


private:
	MOEquivalence & operator = (const MOEquivalence &);

	void	scanString(String);


protected:
	MOEquivalence() {}



INT	n;								// number of equivalences
TEquiv  *equiv;				// list of equivalences

static const EnergyType	degenTolerance;	// tolerance for degeneration
};








#endif
