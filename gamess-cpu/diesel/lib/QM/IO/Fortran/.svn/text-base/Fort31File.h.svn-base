//***********************************************************************
//
//	Name:			Fort31File.h
//
//	Description:	
//
//	Author:			Michael Hanrath
//
//	Version:		0.0
//
//	Date:			10. Jun 1998
//
//***********************************************************************


#ifndef __FORT31FILE_H
#define __FORT31FILE_H

#include "../../../../config.h"

#include "Fort31RecordFormat.h"

#include <stdlib.h>
#include "../../../Container/String.h"

struct Fort31File {
	Fort31File();
	Fort31File(String name);					// auto detect file format
		
	Fort31File(
		String name,
		Fort31RecordFormat	format);
		
	String	name;
	Fort31RecordFormat	format;
};


inline
Fort31File::Fort31File()
{
//	name = NULL;
	format = Fort31RecordFormatUndefined;
}

inline
Fort31File::Fort31File(
	String _name,
	Fort31RecordFormat	_format)
{
	name = _name;
	format = _format;
}		




#endif
