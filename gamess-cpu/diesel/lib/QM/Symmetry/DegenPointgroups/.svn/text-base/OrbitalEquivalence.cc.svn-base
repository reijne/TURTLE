#include "OrbitalEquivalence.h"

#include "../../MO/Iterators/MOEquivalence.h"

#include <algorithm>

using std::find;

OrbitalEquivalence::OrbitalEquivalence(const MOEquivalence &moe)
{
vector<INT>	used;

	if ( &moe )
	{
		for ( INT i=0 ; i<moe.getNumberOfEquivalences() ; ++i )
		{
		vector<INT>	e;
			for ( INT j=0 ; j<moe.getEquiv(i).n ; ++j )
			{
				e.push_back(moe.getEquiv(i).mo[j]);
				used.push_back(moe.getEquiv(i).mo[j]);
			}
			push_back(e);
		}
	}
	for ( INT i=1 ; i<300 ; i++ )
	{
		if ( find(used.begin(), used.end(), i)==used.end() )
		{
		vector<INT>	e;
			e.push_back(i);
			push_back(e);
		}
	}
	
	
}




ConfSet	OrbitalEquivalence::operator () (const FullConf &fc)
{
ConfSet	confSet;
 	const_iterator ii = begin();
	for ( vector<INT>::const_iterator i = fc.begin() ; 
		i!=fc.end() ; ++i , ++ii )
	{
	Occupation	occupation(*i, (*ii).size());
//		cout << occupation << endl;
	ConfSet	cs(occupation, *ii);
//		cout << "::" << cs << endl;
		confSet *= cs;
	}



	return confSet;
}


FullConf OrbitalEquivalence::operator () (const Configuration<MOType> &a)
{
FullConf	fc;
	fc.resize(size());
	for ( unsigned INT i=0 ; i<fc.size() ; ++i )
		fc[i] = 0;
	for ( INT i=0 ; i<a.getNumberOfOpenShells() ; i++ )
	{
	INT	j = getInd(a.getOpenShell(i));
		fc[j]++;
	}
	for ( INT i=0 ; i<a.getNumberOfClosedShells() ; i++ )
	{
	INT	j = getInd(a.getClosedShell(i));
		fc[j] += 2;
	}
	
	return fc;
}

ostream & operator << (ostream &s, const OrbitalEquivalence &f)
{
	for ( OrbitalEquivalence::const_iterator i=f.begin() ; i!=f.end() ; ++i )
	{
		for ( vector<INT>::const_iterator j = i->begin() ; j!=i->end() ; ++j )
			s << *j << " ";
		s << endl;
	}
	return s;
}

