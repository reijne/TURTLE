#include "FullConf.h"


/*
FullConf::FullConf(const Configuration<MOType> &a)
{
}
*/


ostream & operator << (ostream &s, const FullConf &f)
{
	s << "[";
	for ( vector<INT>::const_iterator i=f.begin(); i!=f.end() ; ++i )
		s << *i << " ";
	s << "]";
	return s;
}



