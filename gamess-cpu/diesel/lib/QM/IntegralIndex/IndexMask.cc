//***********************************************************************
//
//	Name:			IndexMask.cc
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




#include "IndexMask.h"




ostream& operator<<(ostream& s, const IndexMask &im)
{
	for ( INT i=0 ; i<im.getNumberOfFlags() ; i++ )
		s << im.getFlag(i) << " ";
	return s;
}
