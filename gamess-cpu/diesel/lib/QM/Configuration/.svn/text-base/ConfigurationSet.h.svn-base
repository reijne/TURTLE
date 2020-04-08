//***********************************************************************
//
//	Name:			ConfigurationSet.h
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

#ifndef __ConfigurationSet_h
#define __ConfigurationSet_h

#include "../../../config.h"

#include <iostream>

#include "../../Container/AVLSet.h"
#include "Configuration.h"

class ConfigurationSet :
	public AVLSet<Configuration<MOType> > {
public:
	ConfigurationSet() { numbering = 0; }
	ConfigurationSet(istream &);
	~ConfigurationSet() {}

	ConfigurationSet(const ConfigurationSet &);
	ConfigurationSet & operator = (const ConfigurationSet &);

	INT operator == (const ConfigurationSet &) const;
	INT operator != (const ConfigurationSet &) const;
	
	void	setNumbering(INT);
	
	void	writeToStream(ostream &) const;

	void projectOnIrrep(INT irrep, const MOIrReps & moirreps);

	friend ostream& operator<<(ostream& s, const ConfigurationSet &);
private:
INT	numbering;
};


inline
INT ConfigurationSet::operator == (const ConfigurationSet &a) const
{
	return ((AVLSet<Configuration<MOType> > *) (ConfigurationSet *) this)
		->operator == ((ConfigurationSet &) a);
}

inline
INT ConfigurationSet::operator != (const ConfigurationSet &a) const
{
	return ((AVLSet<Configuration<MOType> > *) (ConfigurationSet *) this)
		->operator != ((ConfigurationSet &) a);
}


#endif
