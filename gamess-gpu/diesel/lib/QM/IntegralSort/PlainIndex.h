//***********************************************************************
//
//	Name:			PlainIndex.h
//
//	Description:	plain index for integral addressing
//
//	Author:			Michael Hanrath
//
//	Version:		0.0
//
//	Date:			30.08.1996
//
//
//
//
//
//***********************************************************************



#ifndef __PLAININDEX_H
#define __PLAININDEX_H

#include "../../../config.h"


typedef unsigned INT PlainIndexType;

class PlainIndex {
public:
	PlainIndex() {}
	PlainIndex(PlainIndexType ind);
	~PlainIndex() {}

//--------------------------------------------------------------------------

	 PlainIndexType	getIndex() const;

//--------------------------------------------------------------------------

protected:
PlainIndexType	index;
};


inline
PlainIndex::PlainIndex(PlainIndexType ind)
{	index = ind;	}

inline
PlainIndexType	PlainIndex::getIndex() const
{	return	index;	}


inline
unsigned INT	iover2(unsigned INT i)
{	return	i*(i-1) >> 1;	}

#endif
