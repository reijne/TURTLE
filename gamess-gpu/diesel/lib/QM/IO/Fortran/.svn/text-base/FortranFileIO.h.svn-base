//***********************************************************************
//
//	Name:			FortranFileIO.h
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




#ifndef __FORTRANFILEIO_H
#define __FORTRANFILEIO_H

#include "../../../../config.h"


class FortranFileIO {
public:
enum	EndianFormat { Little, Big };
	FortranFileIO(INT n, INT Endian = Little);			// default: intel
	FortranFileIO(const char *s, INT Endian = Little);	// default: intel
	~FortranFileIO();
	
	
	void	rewind();
	void	end();
	void	nextRecord();
	void	previousRecord();
	
	INT	read(void *, INT max);
	INT	read(void *);
	INT	readInRec(void *, INT max);

	INT	write(void *, INT n);
	INT	writeInRec(void *, INT n);
	
	
private:
	void	init(const char *s);
	void	swap(int &);
	
	
char	filename[100];
INT	fd;

INT	bigEndian;

INT	recPos;
};




#endif
