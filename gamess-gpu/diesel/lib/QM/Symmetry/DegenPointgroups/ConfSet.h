#ifndef __ConfSet_h
#define __ConfSet_h

#include <vector>
#include <set>
using std::set;

#include "../../../../config.h"

#include "../../Configuration/Configuration.h"
#include "Occupation.h"

class ConfSet : 
	public set<Configuration<MOType> >
{
public:
//friend OrbitalEquivalence;	

	ConfSet()
	{
	Configuration<MOType>	conf;
		insert(conf);
	}
	
	ConfSet(const Occupation & occ, vector<INT> mos)
	{
		for ( Occupation::const_iterator i = occ.begin() ; i!=occ.end() ; ++i )
		{
		Configuration<MOType>	conf;
			for (
					 vector<INT>::const_iterator j = (*i).begin(),
					 jj = mos.begin()
					 ; j!=(*i).end() ; ++j , ++jj )
				for ( INT k=0 ; k<*j ; ++k )
					conf.create(*jj);
			insert(conf);
		}
	}
	
	ConfSet(const Configuration<MOType> & conf)
	{	insert(conf);	}

	ConfSet & operator +=(ConfSet cs);
	ConfSet & operator *=(ConfSet cs);


	void	projectOnIrRep(const MOIrReps &, IrRep irrep);

	INT getInd(const Configuration<MOType> & conf) const
	{
	INT	ii = 0;
		for ( const_iterator i=begin() ; i!=end() ; ++i , ++ii )
			if ( *i==conf )
				return ii;
	
		return -1;
	}

};

ConfSet operator - (const ConfSet & a, const ConfSet & b);

ostream & operator << (ostream & s, const ConfSet & c);



#endif
