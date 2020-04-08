//***********************************************************************
//
//	Name:			TupelIterator.cc
//
//	Description:	iterators for MO lists
//
//	Author:			Michael Hanrath
//
//	Version:		0.0
//
//	Date:			23.07.1997
//
//
//
//
//
//***********************************************************************




#include "TupelIterator.h"





ostream& operator<<(ostream & s, const TupelIterator & tupel)
{
	s << "[";
	for ( INT i=0 ; i<(tupel.n-tupel.iter)/2 ; i++ )
		s << "2";
	for ( INT i=0 ; i<tupel.iter ; i++ )
		s << "1";
	s << "]";
	return s;
}

