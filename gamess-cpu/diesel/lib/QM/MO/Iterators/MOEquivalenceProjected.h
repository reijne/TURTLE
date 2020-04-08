//***********************************************************************
//
//	Name:			MOEquivalenceProjected.h
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


#ifndef __MOEquivalenceProjected_h
#define __MOEquivalenceProjected_h

#include "../../../../config.h"


#include "../../Configuration/Configuration.h"
#include "MOEquivalence.h"

class ConfigurationSet;
class MRMOs;

class MOEquivalenceProjected : 
	protected MOEquivalence {
public:
	MOEquivalenceProjected(const MOEquivalence &,
		const Configuration<MOType> &,
		const MRMOs *);
	~MOEquivalenceProjected();
	
	
	void	generate(ConfigurationSet &);
	void	generate(ConfigurationSet &,
		 Configuration<MOType> conf,
		 INT level);
	
		static void	distribute(INT n, MOType *mo, INT m);

private:

	struct Iterator {
		Iterator(INT n, INT m);
		~Iterator();
	
	INT	next();

	INT	*ii;			// loop variables
	INT	i;				// loop level
	INT	n;				// # objects
	INT m;				// # containers
	
	};

Configuration<MOType>	conf;
INT	*doubly;
INT	*singly;
IrRep	irrep;
const MRMOs	*mrmos;
};



#endif
