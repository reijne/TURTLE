#ifndef __FullConf_h
#define __FullConf_h

#include <vector>

#include "../../../../config.h"

#include "../../Configuration/Configuration.h"
#include "ConfSet.h"


class FullConf :
	public vector<INT> {
public:
//friend OrbitalEquivalence;	

//	FullConf(const Configuration<MOType> &a);
	
	


	bool operator < (const FullConf & fc) const
	{
		if ( size()<fc.size() )
			return true;
		if ( size()>fc.size() )
			return false;
		for ( unsigned INT i=0 ; i<size() ; ++i )
		{
			if ( (*this)[i]<fc[i] )
				return true;
			if ( (*this)[i]>fc[i] )
				return false;
		}
		return false;
	}


};

ostream & operator << (ostream &s, const FullConf &f);



#endif
