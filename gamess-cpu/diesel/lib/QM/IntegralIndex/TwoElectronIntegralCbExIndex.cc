//***********************************************************************
//
//	Name:			TwoElectronIntegralCbExIndex.cc
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




#include "TwoElectronIntegralCbExIndex.h"




ostream& operator<<(ostream& s, const TwoElectronIntegralCbExIndex & t)
{
	s << "(" << t.getI() << " " << t.getJ() << " | " <<
		t.getK() << " " << t.getL() << "), " << 
		t.getCoulombTriade() << ":" << t.getExchangeTriade();
	return s;
}




