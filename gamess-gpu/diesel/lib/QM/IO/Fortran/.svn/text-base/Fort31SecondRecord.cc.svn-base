//***********************************************************************
//
//	Name:			Fort31SecondRecord.cc
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


#include "Fort31SecondRecord.h"

#include <string>


Fort31SecondRecord::Fort31SecondRecord(
	const char *filename,
	Fort31RecordFormat _format)
{
	fIO = new FortranFileIO(filename);
        format = _format;

	memset(&fort31_2, 0, sizeof(Record));
	
	fIO->nextRecord();
	if ( format==Fort31RecordFormatTRADPT )
		fIO->nextRecord();
		
	fIO->readInRec(NULL , 0);
	fIO->readInRec(fort31_2.nit, 666*sizeof(INT));
	fIO->readInRec(&fort31_2.nn, sizeof(INT));
	fIO->readInRec(fort31_2.lg, 8*sizeof(INT));
	switch ( format ) {
	case Fort31RecordFormatNew:
// 		fIO->readInRec(fort31_2.icf, 6400*sizeof(INT));
		fIO->readInRec(fort31_2.ibal, 8*sizeof(INT));
		fIO->readInRec(fort31_2.itil, 8*sizeof(INT));
		fIO->readInRec(fort31_2.mcomp, 100*sizeof(INT));
		break;
		
	case Fort31RecordFormatTRADPT:
		break;
	}
	
}

Fort31SecondRecord::~Fort31SecondRecord()
{
	delete fIO;
}
