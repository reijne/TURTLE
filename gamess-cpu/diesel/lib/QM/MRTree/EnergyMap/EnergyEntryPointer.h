//***********************************************************************
//
//	Name:			EnergyEntryPointer.h
//
//	Description:	used to garanty the initialization of a pointer
//					with NIL
//
//	Author:			Michael Hanrath
//
//	Version:		0.0
//
//	Date:			10.05.1997
//
//
//
//
//***********************************************************************


#ifndef __EnergyEntryPointer_h
#define __EnergyEntryPointer_h

#include "../../../../config.h"


class EnergyEntry;

struct EnergyEntryPointer {
public:
	EnergyEntryPointer();
	

EnergyEntry	*p;
};


EnergyEntryPointer::EnergyEntryPointer()
{	p = NULL;	}



#endif
