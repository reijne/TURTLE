#include "Occupation.h"



ostream & operator << (ostream & s, const Occupation & occ)
{
	for ( Occupation::const_iterator i = occ.begin() ; i!=occ.end() ; ++i )
	{
		for ( vector<INT>::const_iterator j = (*i).begin() ; j!=(*i).end() ; ++j )
			s << *j << " ";
		s << std::endl;
	}
	return s;
}

Occupation::Occupation(INT max, INT orb)
{
vector<INT>	vv(orb);
INT j = 0;

	for ( INT i=0 ; i<orb ; i++ )
		vv[i] = 0;


	j = orb-1;

	for ( ; ; )
	{
		if ( j==orb-1 )
		{
		INT	sum = 0;
			for ( INT i=0 ; i<orb ; i++ )
				sum += vv[i];

			if ( sum==max )
				push_back(vv);
		}

		if ( vv[j]>=2 )
		{
			vv[j] = 0;
			j--;
			if ( j<0 )
				break;
			continue;
		}

		vv[j]++;
		if ( j<orb-1 )
			j++;
	}
}

