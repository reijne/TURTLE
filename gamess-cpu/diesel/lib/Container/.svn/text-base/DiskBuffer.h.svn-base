//***********************************************************************
//
//	Name:			DiskBuffer.h
//
//	Description:	
//
//	Author:			Michael Hanrath
//
//	Version:		0.0
//
//	Date:			04.12.1996
//
//
//
//
//
//***********************************************************************



#ifndef __DISKBUFFER_H
#define __DISKBUFFER_H

#include "../../config.h"
#include "String.h"


class DiskBuffer {
public:
static const INT	truncateOnOpen = 1;
static const INT	deleteOnClose = 2;
static const INT	noTempDir = 4;
	DiskBuffer(INT ObjectSize, String FileName,
		INT mode = truncateOnOpen | deleteOnClose,
		INT numberOfObjects = 0);
	DiskBuffer(String FileName, INT mode = truncateOnOpen | deleteOnClose);
	~DiskBuffer();


	INT get(INT i, void *) const;
	INT put(INT i, void *);	

	INT	getNumberOfObjects() const;
	INT	getObjectSize() const;

	void	deleteLast();

	void	clear();
	
private:
typedef	INT TObjectSize;
TObjectSize	ObjectSize;
INT	fd;
String	FileName;
INT	mode;
};



inline
INT	DiskBuffer::getObjectSize() const
{
	return	ObjectSize;
}




#endif
