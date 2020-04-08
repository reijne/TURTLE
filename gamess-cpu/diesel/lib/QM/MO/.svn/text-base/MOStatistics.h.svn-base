//***********************************************************************
//
//	Name:			MOStatistics.h
//
//	Description:	
//
//	Author:			Michael Hanrath
//
//	Version:		0.0
//
//	Date:			27.03.1998
//
//
//
//
//
//***********************************************************************

#ifndef __MOStatistics_h
#define __MOStatistics_h

#include "../../../config.h"

#include <iostream>
using std::ostream;

#include "MOType.h"

#include "../RepresentationMatrices/DataTypes.h"
#include "../MRTree/EnergyType.h"

/* FD class ostream; */

class MOStatistics {
public:
	MOStatistics(INT maxMO);
	~MOStatistics();
	
	MOStatistics(const MOStatistics &);


	friend ostream & operator << (ostream &, const MOStatistics &);

private:
	MOStatistics & operator = (const MOStatistics &);

INT	maxMO;					// highest MO number
EnergyType	*e;				// matrix of accumulated MO pair PT energies
INT			*s;				// matrix of number of selected configurations
							// containing certain MO pairs

};





#endif
