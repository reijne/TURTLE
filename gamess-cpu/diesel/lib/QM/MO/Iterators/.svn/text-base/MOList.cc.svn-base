//***********************************************************************
//
//	Name:			MOList.cc
//
//	Description:	MOs classified by internal/external status and
//					irrep to be used as creators or annihilators on
//					internal configuration rests
//
//	Author:			Michael Hanrath
//
//	Version:		0.0
//
//	Date:			05.02.1997
//
//
//
//
//
//***********************************************************************

#include "MOList.h"
#include "../MRMOs.h"

using std::endl;

INT	MOList::getTotalNumber() const
{
INT	number = 1;
	for ( INT i=0 ; i<nTupel ; i++ )
		number *= mrmos->getNumberOfIrRepIntExt(irrep[i], !getIsInternal(i));
	return number;
}




ostream& operator<<(ostream & s, const MOList & molist)
{
	s << "[";
	for ( INT i=0 ; i<molist.nTupel ; i++ )
		s << (i ? ", " : "" ) << "(" << molist.irrep[i] << ", " << 
			(molist.getIsInternal(i) ? "INT" : "ext") << ")";
	s << "]:" << endl;

	return s;
/*
INT	i = 0;
INT	j[molist.nTupel];
	j[i] = molist.mrmos->getIrRepIntExtStart(molist.irrep[i], !molist.isInternal[i]);

	for ( ; ; j[i]++ )
	{
		if ( j[i] >= 
		molist.mrmos->getIrRepIntExtEnd(molist.irrep[i], !molist.isInternal[i]) )
		{
			i--;
			if ( i<0 )
				break;
			else
				continue;
		}
		if ( i<molist.nTupel-1 )
		{
			i++;
			j[i] = molist.mrmos->getIrRepIntExtStart(molist.irrep[i], !molist.isInternal[i]) - 1;
			continue;
		}
		else
		{
			for ( INT k=0 ; k<molist.nTupel ; k++ )
				s << j[k] << " ";
			s << endl;
		}
	}
*/
}



