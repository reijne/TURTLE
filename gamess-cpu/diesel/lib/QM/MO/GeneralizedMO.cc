//***********************************************************************
//
//	Name:			GeneralizedMO.cc
//
//	Description:	
//
//	Author:			Michael Hanrath
//
//	Version:		0.0
//
//	Date:			14.08.1996
//
//
//
//
//
//***********************************************************************




#include "GeneralizedMO.h"



MRMOs * GeneralizedMO::mrmos;


ostream& operator<<(ostream& s, const GeneralizedMO & mo)
{
	if ( mo.getType()==GeneralizedMO::Number )
		s << " " << mo.getMONumber();
	else
	{	s << '(' << mo.getIrRep() << ", ";
		if ( mo.isInternal() )
			s << "INT";
		else
			s << "ext";
		s << ", " << mo.mo << ", " << mo.signature << ')';
	}
	return s;
}
