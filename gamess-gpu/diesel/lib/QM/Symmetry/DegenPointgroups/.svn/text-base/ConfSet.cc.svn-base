#include "ConfSet.h"



ConfSet operator - (const ConfSet & a, const ConfSet & b)
{
ConfSet	diff(a);

	for ( ConfSet::const_iterator i=b.begin() ; i!=b.end() ; ++i )
		diff.erase(*i);
	
	return diff;
}

ostream & operator << (ostream & s, const ConfSet & c)
{
	for ( ConfSet::const_iterator i = c.begin() ; i!=c.end() ; ++i )
		s << *i << endl;
	return s;
}

ConfSet & ConfSet::operator *= (ConfSet cs)
{
ConfSet	c(*this);

	clear();

	for ( ConfSet::const_iterator i = c.begin() ; i!=c.end() ; ++i )
	{
		for ( ConfSet::const_iterator j = cs.begin() ; j!=cs.end() ; ++j )
		{
		Configuration<MOType>	conf(*i);
			conf += *j;
			insert(conf);
		}	
	}	

	return *this;
}

ConfSet & ConfSet::operator += (ConfSet cs)
{
	for ( ConfSet::const_iterator j = cs.begin() ; j!=cs.end() ; ++j )
		insert(*j);
	return *this;
}


void	ConfSet::projectOnIrRep(const MOIrReps &irreps, IrRep irrep)
{
ConfSet	h;
	h.clear();
	for ( const_iterator i = begin() ; i!=end() ; ++i )
		if ( i->calcIrRep(irreps)==irrep )
			h.insert(*i);
	
	*this = h;
}
