//***********************************************************************
//
//	Name:			Fort31FirstRecord.h
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




#ifndef __FORT31FIRSTRECORD_H
#define __FORT31FIRSTRECORD_H

#include "../../../../config.h"

#include "FortranFileIO.h"

#include "Fort31RecordFormat.h"
#include "Fort31File.h"


class Fort31FirstRecord {
public:
	Fort31FirstRecord(Fort31File f31);
	Fort31FirstRecord(const char *filename, Fort31RecordFormat format);
	~Fort31FirstRecord();
        INT getN1el();
	
	INT			getN()			{	return fort31_1.nl;			}
	INT			getNOD()		{	return fort31_1.nodl;		}
	INT			getJOD()		{	return fort31_1.jodl;		}
	INT			getKSUM()		{	return fort31_1.ksuml;		}
	INT			getNBOX()		{	return fort31_1.nboxl;		}
	INT &   	getMJ(INT i)	{	return fort31_1.mjl[i];		}
	INT &   	getKJ(INT i)	{	return fort31_1.kjl[i];		}
	INT &   	getLJ(INT i)	{	return fort31_1.ljl[i];		}
	INT &   	getNJ(INT i)	{	return fort31_1.njl[i];		}
	INT			getNSYM0()		{	return fort31_1.nsym0l;		}
	INT &   	getNTIL(INT i)	{	return fort31_1.ntill[i];	}
	INT &   	getNBAL(INT i)	{	return fort31_1.nball[i];	}
	INT &   	getISYM(INT i)	{	return fort31_1.isyml[i];	}
	INT &   	getJAB(INT i)	{	return fort31_1.jabl[i];	}
	INT			getIORBS()		{	return fort31_1.iorbsl;		}
	INT			getKNU()		{	return fort31_1.knul;		}
	INT &   	getLSYM(INT i)	{	return fort31_1.lsyml[i];	}
	INT &   	getNCOMP(INT i)	{	return fort31_1.ncompl[i];	}
//	double &	getC(INT i)		{	return fort31_1.cl[i];		}
	double		getVNUC()		{	return fort31_1.vnucl;		}
	double		getZERO()		{	return fort31_1.zerol;		}

	
private:
	void	init(
		const char *filename,
		Fort31RecordFormat _format);
		
FortranFileIO	*fIO;
Fort31RecordFormat	format;

struct Record {
	INT			nl;
	INT			nodl;
	INT			jodl;
	INT			ksuml;
	INT			nboxl;
	INT	                mjl[8];
	INT	                kjl[8];
	INT                    	ljl[8];
	INT	                njl[8];
	INT			nsym0l;
	INT             	ntill[8];
	INT	                nball[9];
	INT	                isyml[8];
	INT             	jabl[36];
	INT			iorbsl;
	INT			knul;
	INT             	lsyml[800];
	INT	                ncompl[100];
//	double		cl[200];
	double		vnucl;
	double		zerol;
	} fort31_1;
};




#endif
