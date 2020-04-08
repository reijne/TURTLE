#ifndef __OrbitalEquivalence_h
#define __OrbitalEquivalence_h

#include <vector>

class MOEquivalence;

#include "../../../../config.h"

#include "FullConf.h"
#include "ConfSet.h"

class OrbitalEquivalence :
	public vector<vector<INT> > {
public:
	OrbitalEquivalence() {}
	OrbitalEquivalence(const MOEquivalence &);
	
	
	ConfSet	operator () (const FullConf &);

	FullConf operator () (const Configuration<MOType> &);

private:
	INT getInd(INT mo) const
	{
		for ( unsigned INT i=0 ; i<size() ; ++i )
			for ( unsigned INT j=0 ; j<(*this)[i].size() ; ++j )
				if ( (*this)[i][j]==mo )
					return i;
					
                std::cerr << "mo \"" << mo << "\" not assigned" << endl;
		exit(1);
	}
};


ostream & operator << (ostream &s, const OrbitalEquivalence &f);



#endif
