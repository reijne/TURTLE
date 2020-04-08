#include "../../../config.h"
//***********************************************************************
//
//	Name:			BaseContainer.h
//
//	Description:	
//
//	Author:			Michael Hanrath
//
//	Version:		0.0
//
//	Date:			03.09.1996
//
//
//
//
//
//***********************************************************************

#ifndef __BaseContainer_H
#define __BaseContainer_H


class SharedMemory;

template <class T>
class BaseContainer {
public:
	BaseContainer(INT n);
	~BaseContainer();
	
//-------------------------------------------------------------------------

	INT	getDataSize() const;
	INT	allocate();
	void	check();

	INT	getNumberOfIntegrals() const;	
	INT	getNumberOfLeaves() const;	

	void	setSharedMem(SharedMemory *sharedMem);

//-------------------------------------------------------------------------

protected:
INT	n;
T	**p;
};


#endif
