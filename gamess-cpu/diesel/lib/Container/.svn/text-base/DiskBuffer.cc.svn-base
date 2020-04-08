//***********************************************************************
//
//	Name:			DiskBuffer.cc
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

#include "DiskBuffer.h"
#include "TempDir.h"

#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>
#include <iostream>
#include <string>
#include <unistd.h>
#include <errno.h>
#include <stdlib.h>

using namespace std;

DiskBuffer::DiskBuffer(INT _ObjectSize, String _FileName, INT _mode,
		INT numberOfObjects)
{
	ObjectSize = _ObjectSize;
	mode = _mode;
	
	if ( mode & noTempDir )
		FileName = _FileName;
	else
		FileName = TempDir + "/" + _FileName;
	if ( mode & truncateOnOpen )
		fd = open(FileName.chars(), O_RDWR | O_CREAT | O_TRUNC,
			 S_IRUSR | S_IRGRP | S_IROTH | S_IWUSR);
	else
		fd = open(FileName.chars(), O_RDWR | O_CREAT, 
			S_IRUSR | S_IRGRP | S_IROTH | S_IWUSR);
	write(fd, &ObjectSize, sizeof(TObjectSize));
	if ( numberOfObjects>0 )
	{
	const INT	bufSize = 1<<16;
	char	buf[bufSize];
		memset(buf, 0, bufSize);
		for ( INT i=0 ; i<=ObjectSize*numberOfObjects/bufSize ; i++ )
			if ( i<ObjectSize*numberOfObjects/bufSize )
				write(fd, buf, bufSize);
			else
				write(fd, buf, (ObjectSize*numberOfObjects) % bufSize);
	}

}

DiskBuffer::DiskBuffer(String _FileName, INT _mode)
{
	mode = _mode;
	if ( mode & noTempDir )
		FileName = _FileName;
	else
		FileName = TempDir + "/" + _FileName;
	fd = open(FileName.chars(), O_RDONLY);
	if ( fd!=-1 )
		read(fd, &ObjectSize, sizeof(TObjectSize));
	else
	{
		cerr << "file \"" << _FileName << "\" does not exist." << endl;
		exit(1);
	}
}


DiskBuffer::~DiskBuffer()
{
	close(fd);
	if ( mode & deleteOnClose )
		unlink(FileName);
}



INT	DiskBuffer::get(INT i, void *p) const
{
	lseek(fd, i*ObjectSize + sizeof(TObjectSize), SEEK_SET);
	return (read(fd, p, ObjectSize)!=ObjectSize);
}


INT	DiskBuffer::put(INT i, void *p)
{
LONG_INT	l =	lseek(fd, i*ObjectSize + sizeof(TObjectSize), SEEK_SET);
	if ( (l=write(fd, p, ObjectSize))<ObjectSize )
	{
		cout << "error: no space left on device." << endl;
		cerr << "error: no space left on device." << endl;
		cout << "file: " << FileName << ", errno: " << errno << endl;
		cerr << "file: " << FileName << ", errno: " << errno << endl;
		exit(1);
	}
	return 0;
//	return (write(fd, p, ObjectSize)!=ObjectSize);
}


INT	DiskBuffer::getNumberOfObjects() const
{
struct stat buf;
	fstat(fd, &buf);

	return (buf.st_size - sizeof(TObjectSize))/ObjectSize;
}

void	DiskBuffer::deleteLast()
{
	ftruncate(fd, lseek(fd, 0, SEEK_END) - ObjectSize);
}

void	DiskBuffer::clear()
{
	close(fd);
	fd = open(FileName, O_RDWR | O_CREAT | O_TRUNC);
	write(fd, &ObjectSize, sizeof(TObjectSize));
	chmod(FileName, S_IRUSR | S_IWUSR);
}
