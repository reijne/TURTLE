//***********************************************************************
//
//	Name:			ContainerIterator.h
//
//	Description:	general container iterator
//
//	Author:			Michael Hanrath
//
//	Version:		0.0
//
//	Date:			15.03.1997
//
//
//
//
//
//***********************************************************************


#ifndef __ContainerIterator_h
#define __ContainerIterator_h

#include "../../../../config.h"


#include "../../../Container/Pix.h"



union ContainerIterator {
Pix	pix;
INT	i;

	ContainerIterator()	{}
	ContainerIterator(Pix _pix)	{	pix = _pix;	}
	ContainerIterator(INT _i)	{	i = _i;	}



	INT operator == (const ContainerIterator &);
};

inline
INT ContainerIterator::operator == (const ContainerIterator &a)
{
	return i==a.i;
}

#endif
