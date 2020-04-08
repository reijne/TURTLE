//***********************************************************************
//
//	Name:			TwoElectronIntegralTriadeIndex.cc
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




#include "TwoElectronIntegralTriadeIndex.h"



unsigned short 
	TwoElectronIntegralTriadeIndex::IndexTable[64] =
	//
	//
	//
	//
	{
	0xe464,		// (01|23)
	0xe198,		// (10|23)
	0x0,		// non existent
	0xc912,		// (20|13)
	0x0,		// non existent
	0x0,		// non existent
	0x0,		// non existent
	0x3921,		// (30|12)
	0xd846,		// (02|13)
	0x0,		// non existent
	0xd289,		// (12|03)
	0xc621,		// (21|03)
	0x0,		// non existent
	0x0,		// non existent
	0x0,		// non existent
	0x3612,		// (31|02)
	0x0,		// non existent
	0x0,		// non existent
	0x0,		// non existent
	0x0,		// non existent
	0x0,		// non existent
	0x0,		// non existent
	0x0,		// non existent
	0x0,		// non existent
	0x7889,		// (03|12)
	0x0,		// non existent
	0x7246,		// (13|02)
	0x0,		// non existent
	0x0,		// non existent
	0x0,		// non existent
	0x4e64,		// (23|01)
	0x1e98,		// (32|01)
	0xb498,		// (01|32)
	0xb164,		// (10|32)
	0x0,		// non existent
	0x0,		// non existent
	0x0,		// non existent
	0x8d46,		// (20|31)
	0x0,		// non existent
	0x2d89,		// (30|21)
	0x0,		// non existent
	0x0,		// non existent
	0x0,		// non existent
	0x0,		// non existent
	0x0,		// non existent
	0x0,		// non existent
	0x0,		// non existent
	0x0,		// non existent
	0x9c12,		// (02|31)
	0x0,		// non existent
	0x0,		// non existent
	0x0,		// non existent
	0x9321,		// (12|30)
	0x8789,		// (21|30)
	0x0,		// non existent
	0x2746,		// (31|20)
	0x6c21,		// (03|21)
	0x0,		// non existent
	0x0,		// non existent
	0x0,		// non existent
	0x6312,		// (13|20)
	0x0,		// non existent
	0x4b98,		// (23|10)
	0x1b64,		// (32|10)
	};


ostream& operator<<(ostream& s, const TwoElectronIntegralTriadeIndex & t)
{
	s << "(" << t.getI() << " " << t.getJ() << " | " <<
		t.getK() << " " << t.getL() << "), " << t.getM();
	return s;
}




