//***********************************************************************
//
//	Name:			TypeClass.cc
//
//	Description:	
//
//	Author:			Michael Hanrath
//
//	Version:		0.0
//
//	Date:			27. Jan 1998
//
//***********************************************************************




#include "TypeClass.h"


TypeClass::TypeClass()
{
	n = 0;
}


TypeClass::TypeClass(String _label)
{
	Name = _label;
	n = 0;
}


TypeClass::TypeClass(istream &s)
{
	s >> Name;
	s >> n;
}





String	TypeClass::getName() const
{
	return	Name;
}


void	TypeClass::writeToStream(ostream &s) const
{
	s << Name << " " << n;
}

ostream &	operator << (ostream &s, const TypeClass & a)
{
	s << a.Name;
	return s;
}
