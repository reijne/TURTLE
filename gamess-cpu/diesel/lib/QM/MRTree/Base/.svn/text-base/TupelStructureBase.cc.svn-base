//***********************************************************************
//
//	Name:			TupelStructureBase.cc
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

#include "TupelStructureBase.h"
#include "extMOsBase.h"
#include "InternalConfsBase.h"
#include "NExternalsBase.h"

#include <stdio.h>


TupelStructureBase::TupelStructureBase(istream &s) :
		Configuration<MOType>(s)
{
//	cout << "TupelStructureBase(istream &s)" << endl;
	s >> irrep;
	s >> ReferenceFlag;

//	cout << "(((((((((((((((((((" << endl;
//	cout << *this;
//	cout << ")))))))))))))))))))" << endl;
}



TupelStructureBase::TupelStructureBase(IrRep _irrep, 
		INT NumberOfOpenShells, INT NumberOfClosedShells, MOType *pShells,
		INT _ReferenceFlag) :
	Configuration<MOType>(
		NumberOfOpenShells,
		NumberOfClosedShells,
		pShells)
{
	irrep = _irrep;
	ReferenceFlag = _ReferenceFlag;
}


TupelStructureBase::TupelStructureBase(IrRep _irrep, 
		Configuration<MOType> conf,
		INT _ReferenceFlag) :
	Configuration<MOType>(conf)
{
	irrep = _irrep;
	ReferenceFlag = _ReferenceFlag;
}


TupelStructureBase::~TupelStructureBase()
{
}


/*
void	TupelStructureBase:setNumber(INT n)
{
	if ( number == n )
		return;
	delete moContainer;
	number = n;
}
*/

void	TupelStructureBase::writeToStream(ostream & s) const
{
	Configuration<MOType>::writeToStream(s);
	s << endl;
	s << irrep << endl;
	s << ReferenceFlag << endl;
}

ostream& operator<<(ostream & s, const TupelStructureBase &base)
{
	s << ((Configuration<MOType> &) base) << endl;
	s << base.irrep << endl;
	s << base.ReferenceFlag << endl;
	return s;	
}


istream& operator>>(istream & s, TupelStructureBase &base)
{
	s >> ((Configuration<MOType> &) base);
	s >> base.irrep;
	s >> base.ReferenceFlag;
	return s;	
}

