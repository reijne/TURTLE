//***********************************************************************
//
//	Name:			OrderedSetElement.h
//
//	Description:	abstract class containing operators supported by
//					elements of an ordered set
//
//	Author:			Michael Hanrath
//
//	Version:		0.0
//
//	Date:			02.06.1997
//
//
//
//
//
//***********************************************************************


#ifndef __OrderedSetElement_h
#define __OrderedSetElement_h

#include "../../../../config.h"



class OrderedSetElement {
public:

	//--------------------------------------------------------------------------

	virtual INT operator == (const OrderedSetElement &) = 0;
	virtual INT operator != (const OrderedSetElement &) = 0;
	virtual INT operator <= (const OrderedSetElement &) = 0;

	//----------------------------------------------------------------




};



#endif
