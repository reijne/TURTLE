//***********************************************************************
//
//	Name:			SharedMemory.h
//
//	Description:	
//
//	Author:			Michael Hanrath
//
//	Version:		0.0
//
//	Date:			20.05.1997
//
//
//
//
//
//***********************************************************************

#ifndef __SharedMemory_h
#define __SharedMemory_h

#include "../../../config.h"

#include <stdio.h>

typedef INT key_T;

class SharedMemory {
public:
enum Mode { StandAlone, Master, Slave };

	SharedMemory(INT size, Mode mode);		// constructor for 
											// StandAlone/Master mode
											
	SharedMemory(key_T key, INT size);		// constructor for Slave mode

	~SharedMemory();


	void	setSlave();
	
	
	void *	allocate(INT n, INT size);
	void	deallocate(void *);

	key_T	getKey() const;

private:
	void	freeShm();
	
char	*p;							// pointer to allocated memory
INT	size;							// size of memory block
Mode	mode;						// flag if master or slave
key_T	key;						// shared memory system key
INT	id;								// shared memory system id
INT	tempSize;						// temporary size

static key_T	nextKey;			// next key to use
};



#endif
