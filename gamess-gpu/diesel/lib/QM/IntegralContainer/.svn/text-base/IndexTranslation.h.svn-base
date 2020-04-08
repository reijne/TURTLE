#include "../../../config.h"
//***********************************************************************
//
//	Name:			IndexTranslation.h
//
//	Description:	
//
//	Author:			Michael Hanrath
//
//	Version:		0.0
//
//	Date:			03.09.1996
//
//
//
//
//
//***********************************************************************

#ifndef __INDEXTRANSLATION_H
#define __INDEXTRANSLATION_H


#include "../MO/MOType.h"


class IndexTranslation {
public:
	IndexTranslation(INT maxMO);
	~IndexTranslation();
	
//------------------------------------------------------------------------

	const INT *	getIndex(INT i) const;
	
	INT	isPointer(const INT * p) const;

//------------------------------------------------------------------------

private:
INT	*Index[5];		//	0: 0
					//	1: i
					//	2: i*(i+1)/2
					//	3: i*(i+1)*(i+2)/6
					//	4: i*(i+1)*(i+2)*(i+3)/24
};



inline
const INT *	IndexTranslation::getIndex(INT i) const
{	return Index[i];	}


inline
INT	IndexTranslation::isPointer(const INT * p) const
{
	for ( INT i=0 ; i<5 ; i++ )
		if ( p == Index[i] )
			return 1;
	return 0;
}


#endif
