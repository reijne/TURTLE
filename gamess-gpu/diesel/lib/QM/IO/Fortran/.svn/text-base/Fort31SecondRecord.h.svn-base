//***********************************************************************
//
//	Name:			Fort31SecondRecord.h
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




#ifndef __FORT31SECONDRECORD_H
#define __FORT31SECONDRECORD_H

#include "../../../../config.h"

#include "FortranFileIO.h"

#include "Fort31RecordFormat.h"


class Fort31SecondRecord {
public:
	Fort31SecondRecord(const char *filename, Fort31RecordFormat format);
	~Fort31SecondRecord();
	
	
	INT	&		getNIT(INT i)	{	return fort31_2.nit[i];		}
	INT			getNN()			{	return fort31_2.nn;			}
	INT &	                getLG(INT i)	{	return fort31_2.lg[i];		}
//	INT &	                getICF(INT i)	{	return fort31_2.icf[i];		}
	INT &	                getIBAL(INT i)	{	return fort31_2.ibal[i];	}
	INT &	                getITIL(INT i)	{	return fort31_2.itil[i];	}
	INT &	                getMCOMP(INT i)	{	return fort31_2.mcomp[i];	}

	
	
private:
FortranFileIO	*fIO;
Fort31RecordFormat	format;
struct Record {
	INT			nit[666];
	INT			nn;
	INT             	lg[8];
//	INT             	icf[6400];
	INT             	ibal[8];
	INT             	itil[8];
	INT	                mcomp[100];
	} fort31_2;
};




#endif
