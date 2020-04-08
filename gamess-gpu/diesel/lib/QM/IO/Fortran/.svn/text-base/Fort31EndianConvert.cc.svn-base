//***********************************************************************
//
//	Name:			Fort31EndianConvert.cc
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


#include "Fort31EndianConvert.h"

#include <stdio.h>




Fort31EndianConvert::Fort31EndianConvert(
		const char *inFileName,
		const char *outFileName,
		Direction _dir,
		Fort31RecordFormat _format)
{
	dir = _dir;
	format = _format;

INT	*buf = new INT[1000000];



FortranFileIO	*fin = new FortranFileIO(inFileName, dir==BigToLittle);
FortranFileIO	*fout = new FortranFileIO(outFileName, dir==LittleToBig);

INT	l = 0;
	switch ( format )
	{
	case Fort31RecordFormatOld:
	{

	const INT	FirstEndian = 7;
	TData	FirstRecordEndian[FirstEndian] = {
		{ 5, 4 },
		{ 32, 2 },
		{ 1, 4 },
		{ 61, 2 },
		{ 2, 4 },
		{ 900, 2 },
		{ 2677, 8 }
		};	

	const INT	SecondEndian = 4;
	TData	SecondRecordEndian[SecondEndian] = {
		{ 667, 4 },
		{ 8, 2 },
		{ 6400, 8 },
		{ 6516, 2 }
		};	


	const INT	DataEndian = 1;
	TData	DataRecordEndian[DataEndian] = {
		{ 2000, 4 }
		};	


	const INT	LastEndian = 3;
	TData	LastRecordEndian[LastEndian] = {
		{ 5000, 8 },
		{ 10000, 4 },
		{ 1, 8 }
		};	


		fin->read(buf);
		fout->write(buf, swap((char*) buf, FirstRecordEndian, FirstEndian));

		fin->read(buf);
		if ( dir==LittleToBig )
			l = buf[666];
		fout->write(buf, swap((char*) buf, SecondRecordEndian, SecondEndian));
		if ( dir==BigToLittle )
			l = buf[666];
		printf("%d integrals\n", l);

		for ( INT i=0 ; i<(l+1999)/2000 ; i++ )
		{	fin->read(buf);
			fout->write(buf, swap((char*) buf, DataRecordEndian, DataEndian));
		}

		fin->read(buf);
	fout->write(buf, swap((char*) buf, LastRecordEndian, LastEndian));
	}
		break;


	case Fort31RecordFormatNew:
	{

	const INT	FirstEndian = 7;
	TData	FirstRecordEndian[FirstEndian] = {
		{ 5, 4 },
		{ 32, 2 },
		{ 1, 4 },
		{ 61, 2 },
		{ 2, 4 },
		{ 900, 2 },
		{ 202, 8 }
		};	

	const INT	SecondEndian = 3;
	TData	SecondRecordEndian[SecondEndian] = {
		{ 667, 4 },
		{ 8, 2 },
		{ 6516, 2 }
		};	

		fin->read(buf);
		fout->write(buf, swap((char*) buf, FirstRecordEndian, FirstEndian));

		fin->read(buf);
		if ( dir==LittleToBig )
			l = buf[666];
		fout->write(buf, swap((char*) buf, SecondRecordEndian, SecondEndian));
		if ( dir==BigToLittle )
			l = buf[666];
		printf("%d integrals\n", l);


	const INT	DataEndian = 1;
	TData	DataRecordEndian[DataEndian] = {
		{ 2000, 4 }
		};	


		for ( INT i=0 ; i<(l+1999)/2000 ; i++ )
		{	fin->read(buf);
			fout->write(buf, swap((char*) buf, DataRecordEndian, DataEndian));
		}

	const INT	LastEndian1 = 1;
	TData	LastRecordEndian1[LastEndian1] = {
		{ 31376, 8 }
		};	
		fin->read(buf);
		fout->write(buf, swap((char*) buf, LastRecordEndian1, LastEndian1));

	const INT	LastEndian2 = 1;
	TData	LastRecordEndian2[LastEndian2] = {
		{ 31376, 4 }
		};	
		fin->read(buf);
		fout->write(buf, swap((char*) buf, LastRecordEndian2, LastEndian2));

	const INT	LastEndian3 = 1;
	TData	LastRecordEndian3[LastEndian3] = {
		{ 31376, 4 }
		};	
		fin->read(buf);
		fout->write(buf, swap((char*) buf, LastRecordEndian3, LastEndian3));

	const INT	LastEndian4 = 1;
	TData	LastRecordEndian4[LastEndian4] = {
		{ 1, 8 }
		};	
		fin->read(buf);
		fout->write(buf, swap((char*) buf, LastRecordEndian4, LastEndian4));
	}
		break;
		

	case Fort31RecordFormatTRADPT:
	{

	const INT	FirstEndian = 6;
	TData	FirstRecordEndian[FirstEndian] = {
		{ 5, 4 },
		{ 32, 2 },
		{ 1, 4 },
//		{ 68, 2 },
		{ 61, 2 },
		{ 1, 4 },
		{ 2, 8 }
		};	

	const INT	SecondEndian = 2;
	TData	SecondRecordEndian[SecondEndian] = {
		{ 667, 4 },
		{ 8, 2 }
		};	


	const INT	DataEndian = 1;
	TData	DataRecordEndian[DataEndian] = {
		{ 2000, 4 }
		};	

	TData	Data2RecordEndian[DataEndian] = {
		{ 2000, 8 }
		};	

	TData	Data3RecordEndian[DataEndian] = {
		{ 4000, 4 }
		};	


	const INT	LastEndian = 1;
	TData	LastRecordEndian[LastEndian] = {
		{ 1, 8 }
		};	

		fout->write(buf, fin->read(buf));

		fin->read(buf);

	INT nboxl = 0;
		if ( dir==LittleToBig )
			nboxl = buf[4];
		fout->write(buf, swap((char*) buf, FirstRecordEndian, FirstEndian));
		if ( dir==BigToLittle )
			nboxl = buf[4];

		fin->read(buf);
		if ( dir==LittleToBig )
			l = buf[666];
		fout->write(buf, swap((char*) buf, SecondRecordEndian, SecondEndian));
		if ( dir==BigToLittle )
			l = buf[666];
		printf("%d integrals\n", l);

		for ( INT i=0 ; i<(l+1999)/2000 ; i++ )
		{	fin->read(buf);
			fout->write(buf, swap((char*) buf, DataRecordEndian, DataEndian));
		}

		l = nboxl*(nboxl+1)/2;
		for ( INT i=0 ; i<l/2000 ; i++ )
		{	fin->read(buf);
			fout->write(buf, swap((char*) buf, Data2RecordEndian, DataEndian));
		}
		if ( l%2000 )
		{
			fin->read(buf);
			swap((char*) buf, Data2RecordEndian, DataEndian);
			fout->write(buf, (l % 2000)*8);
		}

		for ( INT i=0 ; i<l/2000 ; i++ )
		{	fin->read(buf);
			fout->write(buf, swap((char*) buf, Data3RecordEndian, DataEndian));
		}
		if ( l%2000 )
		{
			fin->read(buf);
			swap((char*) buf, Data3RecordEndian, DataEndian);
			fout->write(buf, (l % 2000)*8);
		}




		printf("%d\n", fin->read(buf));
		fout->write(buf, swap((char*) buf, LastRecordEndian, LastEndian));
	}
		break;
		
	}


}


INT	Fort31EndianConvert::swap(char *p, TData *data, INT n)
{
INT	l = 0;
//char   *p = (char *) p1;        

	for ( INT i=0 ; i<n ; i++ )
	{	l += data[i].n*data[i].size;
		for ( INT j=0 ; j<data[i].n ; j++ )
		{	for ( INT k=0 ; k<data[i].size/2 ; k++ )
			{	
			char	h = p[k];
				p[k] = p[data[i].size-1-k];
				p[data[i].size-1-k] = h;
				
			}
			p += data[i].size;
		}
	}
//	printf("recsize: %d\n", l);
	return l;
}



