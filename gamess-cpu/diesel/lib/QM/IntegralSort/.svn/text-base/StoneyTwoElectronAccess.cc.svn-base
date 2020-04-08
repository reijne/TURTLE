//***********************************************************************
//
//	Name:			StoneyTwoElectronAccess.cc
//
//	Description:	implements
//						1.	direct access to stoney integrals
//						2.	
//							(for sorting of integrals)
//
//	Author:			Michael Hanrath
//
//	Version:		0.0
//
//	Date:			30.08.1996
//
//
//
//
//
//***********************************************************************

#include "StoneyTwoElectronAccess.h"

#include "../IO/Fortran/FortranFileIO.h"

#include <stdlib.h>
#include <iostream>

using namespace std;

StoneyTwoElectronAccess::StoneyTwoElectronAccess(
	const char *_Fort31FileName, INT _maxListEntries)
{
	Fort31FileName = _Fort31FileName;
	maxListEntries = _maxListEntries;
	list = new TList[maxListEntries];
	n = 0;
}


StoneyTwoElectronAccess::~StoneyTwoElectronAccess()
{
	delete list;
}



StoneyTwoElectronType	StoneyTwoElectronAccess::getIntegral(
	const PlainIndex &) const
{
}


static int	TListCmp(PlainIndexType *p1, PlainIndexType *p2)
{
	if ( *p1==*p2 )
		return 0;
	if ( *p1>*p2 )
		return 1;
	return -1;
}


void	StoneyTwoElectronAccess::getRequestedIntegrals(StoneyTwoElectronType *p)
{
	qsort(list, n, sizeof(TList), TListCmp);

//	for ( INT i=0 ; i<n ; i++ )
//		cout << i << ". " << list[i].request << "  " << list[i].memIndex << endl;

FortranFileIO	fort31(Fort31FileName);

	fort31.nextRecord();
	fort31.nextRecord();

INT	alt, neu;
const INT	recSize = 2000;
StoneyTwoElectronType	buf[recSize];

	alt = -1;
	for ( INT i=0 ; i<n ; i++ )
	{	neu = list[i].request / recSize;
		while ( neu!=alt )
		{	fort31.read(buf);
//			for ( INT j=0 ; j<10 ; j++ ) cout << buf[j] << " ";
			alt++;
		}
		p[list[i].memIndex] = buf[list[i].request - alt*recSize];
	}
	

}


