//***********************************************************************
//
//	Name:			ConfigurationSet.cc
//
//	Description:	
//
//	Author:			Michael Hanrath
//
//	Version:		0.0
//
//	Date:			10.06.1997
//
//
//
//
//
//***********************************************************************

#include "ConfigurationSet.h"

#include "../../Container/String.h"
#include <sstream>
using std::stringstream;

#include <iomanip>


ConfigurationSet::ConfigurationSet(istream &s)
{
	while ( s )
	{
	Configuration<MOType>	p(s);
		if ( p.getNumberOfOpenShells() && p.getNumberOfClosedShells() )
			add(p);
	}
	numbering = 0;
}

ConfigurationSet::ConfigurationSet(const ConfigurationSet &set)
{
Pix	i = set.first();
	while(i)
	{
	Configuration<MOType>	n(set(i));
		add(n);
		set.next(i);
	}
	numbering = set.numbering;
}

ConfigurationSet & ConfigurationSet::operator = (const ConfigurationSet &set)
{
	clear();
Pix	i = set.first();
	while(i)
	{
	Configuration<MOType>	n(set(i));
		add(n);
		set.next(i);
	}
	numbering = set.numbering;
	return *this;
}


void ConfigurationSet::projectOnIrrep(INT irrep, const MOIrReps & moirreps)
{
ConfigurationSet	h;
Pix	i = first();
	while ( i )
	{
		if ( (*this)(i).calcIrRep(moirreps) == irrep )
			h.add((*this)(i));
		next(i);
	}
	*this = h;
}


void	ConfigurationSet::writeToStream(ostream &s) const
{
Pix	i = first();
INT	j = 0;
String	l;
INT maxl = 0;
	while ( i )
	{
	stringstream	str(stringstream::in | stringstream::out);
	
		str << (*this)(i) << '\0';
		
		if ( str.tellp()>maxl )
			maxl = str.tellp();
		next(i);
	}

	i = first();
	while ( i )
	{
	stringstream	str(stringstream::in | stringstream::out);
	
		str << (*this)(i) << '\0';
		s << str.str() << std::setw(maxl-str.tellp()) << " ";
		
		if ( numbering )
			s << "\t\t# " << ++j;
		s << endl;	
		next(i);
	}
}



ostream& operator<<(ostream& s, const ConfigurationSet &set)
{
	set.writeToStream(s);
	return s;
}


void	ConfigurationSet::setNumbering(INT f)
{
	numbering = f;
}




