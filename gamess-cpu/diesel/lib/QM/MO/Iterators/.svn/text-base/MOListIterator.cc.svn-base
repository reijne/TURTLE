//***********************************************************************
//
//	Name:			MOListIterator.cc
//
//	Description:	iterators for MO lists
//
//	Author:			Michael Hanrath
//
//	Version:		0.0
//
//	Date:			26.03.1997
//
//
//
//
//
//***********************************************************************




#include "MOListIterator.h"



ostream& operator<<(ostream & s, const MOListIterator & moiter)
{
	s << "mode: ";
	switch ( moiter.mode ) {
	case MOListIterator::noIndex:
		s << "noIndex" << endl;
		break;

	case MOListIterator::oneIndex:
		s << "oneIndex" << ", ";
		s << "i=" << moiter.i << " in [" << moiter.si << "..." << moiter.ei << "]";
		break;
		
	case MOListIterator::twoIndex:
		s << "twoIndex" << ", ";
		s << "i=" << moiter.i << " in [" << moiter.si << "..." << moiter.ei << "]";
		s << ", j=" << moiter.j << " in [" << moiter.sj << "..." << moiter.ej << "]";
		break;
		
	case MOListIterator::lowerTriangular:
		s << "lowerTriangular" << ", ";
		s << "i=" << moiter.i <<  ", j=" << moiter.j 
			<< " in [" << moiter.si << "..." << moiter.ej << "]";
		s << ", i<j";
		break;
	}
	s << ", ";
	
	s << "N=" << moiter.getN();
		
	return s;		
}
