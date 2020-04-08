//***********************************************************************
//
//	Name:			MRTreeIterator.h
//
//	Description:	iterator for MR-tree
//
//	Author:			Michael Hanrath
//
//	Version:		0.0
//
//	Date:			17.03.1997
//
//
//
//
//
//***********************************************************************


#ifndef __MRTreeIterator_h
#define __MRTreeIterator_h

#include "../../../config.h"


#include "Container/ContainerIterator.h"


class	extMOsBase;


struct MRTreeIterator {

	
	INT operator == (const MRTreeIterator &);


	friend ostream& operator<<(ostream & s, const MRTreeIterator &);


ContainerIterator	Citer[4];
extMOsBase	*extMOs;
INT	n;							// counts backward
INT	i;							// counts forward
};



inline
INT MRTreeIterator::operator == (const MRTreeIterator &a)
{
	return 
		Citer[0]==a.Citer[0] &&
		Citer[1]==a.Citer[1] &&
		Citer[2]==a.Citer[2] &&
		Citer[3]==a.Citer[3];
}

#endif
