//***********************************************************************
//
//	Name:			InternalConfsBase.cc
//
//	Description:	
//
//	Author:			Michael Hanrath
//
//	Version:		0.0
//
//	Date:			09.03.1997
//
//
//
//
//
//***********************************************************************

#include "InternalConfsBase.h"
#include "TupelStructureBase.h"

#include <stdio.h>

InternalConfsBase::InternalConfsBase(INT _nExt)
{
	nExt = _nExt;
}

InternalConfsBase::InternalConfsBase(istream &s)
{
	s >> nExt;
//	cout << "InternalConfsBase::InternalConfsBase(istream &s):" << endl;
//	cout << "nExt=" << nExt << endl;
}

/*
InternalConfsBase::InternalConfsBase()
{
}


InternalConfsBase::~InternalConfsBase()
{
}
*/

/*
void	InternalConfs::setNumberOfConfigurations(INT n)
{
	number = n;
	mainConfiguration = new TupelStructure[number];
}
*/

void	InternalConfsBase::writeToStream(ostream & s) const
{
	s << nExt << endl;
}


ostream& operator<<(ostream & s, const InternalConfsBase &base)
{
	s << base.nExt << endl;
	return s;	
}

/*
istream& operator>>(istream & s, InternalConfsBase &base)
{
	return s;	
}
*/
