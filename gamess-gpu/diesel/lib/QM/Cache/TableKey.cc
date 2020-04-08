//***********************************************************************
//
//	Name:	TableKey.cc
//
//	Description:	key for addressing representation matrices
//
//	Author:	Michael Hanrath
//
//	Version:	0.0
//
//	Date:	21.08.1996
//
//
//
//
//***********************************************************************




#include "TableKey.h"


ostream& operator<<(ostream& s, const TableKey & key)
{
	for ( INT i=31 ; i>=0 ; i-- )
	{	if ( key.getKey() & (1<<i) )
			s << "1";
		else
			s << "0";
		if ( i==29 || i==26 || i==24 || i==20 || i==10 )
			s << " ";
	}
	return s;
}




