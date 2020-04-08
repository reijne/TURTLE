//***********************************************************************
//
//	Name:			PackedIntStorage.cc
//
//	Description:	
//
//	Author:			Michael Hanrath
//
//	Version:		0.0
//
//	Date:			19.10.1998
//
//
//
//
//
//***********************************************************************




#include "PackedIntStorage.h"

#include <math.h>

#include <iostream>


PackedIntStorage::PackedIntStorage(INT n, INT _ldBits)
{
	ldBits = _ldBits;
INT	bits = 1 << ldBits;
INT	bitsPerInt = sizeof(INT)*8;
INT	ldBitsPerInt = (INT) floor(log((double)bitsPerInt)/log((double)2) + 0.5);

	mask = (1 << bits) - 1;


	binDiv = ldBitsPerInt - ldBits;
	
	
	binMod = (1 << binDiv) - 1;

	p = new unsigned INT[bits * n / bitsPerInt+1];
}



PackedIntStorage::~PackedIntStorage()
{
	delete p;
}
