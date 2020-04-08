#include "FullConfC2.h"
#include <iomanip>



ostream & operator << (ostream &s, const FullConfC2 &f)
{
	s << "{" << endl;
	s << "[";
	for ( vector<INT>::const_iterator i=f.begin(); i!=f.end() ; ++i )
		s << *i << " ";
	s << "]" << endl;
	
	for ( vector<vector<double> >::const_iterator i=f.c2.begin(); i!=f.c2.end() ; ++i )
	{
		for ( vector<double>::const_iterator j=i->begin(); j!=i->end() ; ++j )
			s << std::setw(12) << *j << " ";
		s << endl;
	}
	
	s << "}" << endl;
	return s;
}



