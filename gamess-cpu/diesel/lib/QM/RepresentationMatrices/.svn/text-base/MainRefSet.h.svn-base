//***********************************************************************
//
//	Name:			MainRefSet.h
//
//	Description:	
//
//	Author:			Michael Hanrath
//
//	Version:		0.0
//
//	Date:			25.08.1999
//
//
//
//
//
//***********************************************************************

#ifndef __MainRefSet_h
#define __MainRefSet_h

#include "../../../config.h"
#include <vector>
using std::vector;

#include "../Configuration/ConfigurationSet.h"

class MainRefSet {
public:
	MainRefSet(INT roots, double limit) : roots(roots), limit(limit) {}

	void	add(Configuration<MOType>, const double *ci2);
	
	void	calc();
	
	ConfigurationSet	getMainRefSet() const;

struct  Data {
        Configuration<MOType>   conf;
        vector<double>  ci2;
        };


private:
INT	roots;
double	limit;
ConfigurationSet	confs;


vector<Data>	data;
};

inline
ConfigurationSet	MainRefSet::getMainRefSet() const
{	return	confs;	}



#endif
