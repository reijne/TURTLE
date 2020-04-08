//***********************************************************************
//
//	Name:			TAccessList.h
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




#ifndef __TACCESSLIST_H
#define __TACCESSLIST_H

#include "../../../config.h"


struct TAccessList {
	void	*key;			// key
	unsigned INT	used;	// time code of last access
	unsigned INT	size;	// size of linked object in byte
	void	*p;			// pointer to object array
	};


#endif
