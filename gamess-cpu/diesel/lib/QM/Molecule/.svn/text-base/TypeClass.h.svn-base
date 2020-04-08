//***********************************************************************
//
//	Name:			TypeClass.h
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

#ifndef _TYPECLASS_H
#define _TYPECLASS_H

#include "../../../config.h"

#include "../../Container/String.h"

#include <iostream>
using std::ostream;
using std::istream;

class TypeClass {
public:
	TypeClass();
	TypeClass(String label);
	TypeClass(istream &s);
	
	INT	operator == (const TypeClass & a) const;
	INT	operator != (const TypeClass & a) const;
	
	void	writeToStream(ostream &s) const;
	
	String	getName() const;

	friend ostream &	operator << (ostream &s, const TypeClass &);

protected:
String	Name;
INT	n;
};


inline
INT	TypeClass::operator == (const TypeClass & a) const
{	return Name==a.Name;	}

inline
INT	TypeClass::operator != (const TypeClass & a) const
{	return Name!=a.Name;	}

#endif
