//***********************************************************************
//
//	Name:			Fort31FirstRecord.cc
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

#include "Fort31FirstRecord.h"

#include <string>
#include <iostream>

Fort31FirstRecord::Fort31FirstRecord(Fort31File f31)
{
	init(f31.name, f31.format);
}

Fort31FirstRecord::Fort31FirstRecord(
	const char *filename,
	Fort31RecordFormat _format)
{
	init(filename, _format);
}

void	Fort31FirstRecord::init(
	const char *filename,
	Fort31RecordFormat _format)
{
	fIO = new FortranFileIO(filename);
	format = _format;

	
	if ( format==Fort31RecordFormatTRADPT )
		fIO->nextRecord();
		
	memset(&fort31_1, 0, sizeof(Record));

	fIO->readInRec(NULL , 0);
	fIO->readInRec(&fort31_1.nl, sizeof(INT));
	fIO->readInRec(&fort31_1.nodl, sizeof(INT));
	fIO->readInRec(&fort31_1.jodl, sizeof(INT));
	fIO->readInRec(&fort31_1.ksuml, sizeof(INT));
	fIO->readInRec(&fort31_1.nboxl, sizeof(INT));
	fIO->readInRec(fort31_1.mjl, 8*sizeof(INT));
	fIO->readInRec(fort31_1.kjl, 8*sizeof(INT));
	fIO->readInRec(fort31_1.ljl, 8*sizeof(INT));
	fIO->readInRec(fort31_1.njl, 8*sizeof(INT));
	fIO->readInRec(&fort31_1.nsym0l, sizeof(INT));
	fIO->readInRec(fort31_1.ntill, 8*sizeof(INT));
	fIO->readInRec(fort31_1.nball, 9*sizeof(INT));
	fIO->readInRec(fort31_1.isyml, 8*sizeof(INT));
	fIO->readInRec(fort31_1.jabl, 36*sizeof(INT));
	fIO->readInRec(&fort31_1.iorbsl, sizeof(INT));
        switch ( format ) {
	case Fort31RecordFormatNew:
		fIO->readInRec(&fort31_1.knul, sizeof(INT));
		fIO->readInRec(fort31_1.lsyml, 800*sizeof(INT));
		fIO->readInRec(fort31_1.ncompl, 100*sizeof(INT));
//		fIO->readInRec(fort31_1.cl, 200*sizeof(double));
		fIO->readInRec(&fort31_1.vnucl, sizeof(double));
		fIO->readInRec(&fort31_1.zerol, sizeof(double));
		break;
		
	case Fort31RecordFormatTRADPT:
		fIO->readInRec(&fort31_1.vnucl, sizeof(double));
		fIO->readInRec(&fort31_1.zerol, sizeof(double));
		break;
		
	}
	

}

Fort31FirstRecord::~Fort31FirstRecord()
{
	delete fIO;
}

INT Fort31FirstRecord::getN1el()
{
   INT i;
   INT n1el = 0;
   for (i=0;i<8;i++) {
     n1el = n1el + fort31_1.ljl[i];
   }
   n1el=n1el*(n1el+1)/2;
   return n1el;
}
