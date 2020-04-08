//***********************************************************************
//
//	Name:			Fort31EndianConvert.h
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




#ifndef __FORT31ENDIANCONVERT_H
#define __FORT31ENDIANCONVERT_H

#include "../../../../config.h"



#include "FortranFileIO.h"

#include "Fort31FirstRecord.h"

class Fort31EndianConvert {
public:

enum	Direction { LittleToBig, BigToLittle };
	Fort31EndianConvert(
		const char *inFileName,
		const char *outFileName,
		Direction dir,
		Fort31RecordFormat format);
	~Fort31EndianConvert() {}
	
	
private:
struct	TData {
	INT	n;
	INT size;
	};

	INT	swap(char  *, TData *, INT);

Direction dir;
Fort31RecordFormat format;
};



#endif
