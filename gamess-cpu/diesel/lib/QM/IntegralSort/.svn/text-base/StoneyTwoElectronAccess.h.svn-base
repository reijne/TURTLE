//***********************************************************************
//
//	Name:			StoneyTwoElectronAccess.h
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



#ifndef __STONEYTWOELECTRONACCESS_H
#define __STONEYTWOELECTRONACCESS_H

#include "../../../config.h"



#include "PlainIndex.h"

typedef double StoneyTwoElectronType;

class StoneyTwoElectronAccess {
public:
	StoneyTwoElectronAccess(const char *Fort31FileName, INT maxListEntries);
	~StoneyTwoElectronAccess();
	
//--------------------------------------------------------------------------

	INT	getNumberOfMaxListEntries() const;
	INT	getNumberOfListEntries() const;
	
//--------------------------------------------------------------------------

//	1. direct access:

	//	immedeate access to an integral (very inefficient)
	StoneyTwoElectronType	getIntegral(const PlainIndex &) const;

//--------------------------------------------------------------------------

//	2.	access through entry list
	
	//	put request for integral in list
	void	putInList(const PlainIndex &);

	//	put the requested integrals in the desired ordering in p
	void	getRequestedIntegrals(StoneyTwoElectronType *p);

	void	clearList();
	
//--------------------------------------------------------------------------
	
private:
const char	*Fort31FileName;
INT	maxListEntries;
INT	n;
struct TList {
	PlainIndexType	request;
	PlainIndexType	memIndex;
	};
TList	*list;
};



inline
INT	StoneyTwoElectronAccess::getNumberOfMaxListEntries() const
{	return	maxListEntries;	}

inline
INT	StoneyTwoElectronAccess::getNumberOfListEntries() const
{	return	n;	}

inline
void	StoneyTwoElectronAccess::putInList(const PlainIndex & plainIndex)
{
	list[n].request = plainIndex.getIndex();
	list[n].memIndex = n++;
}

inline
void	StoneyTwoElectronAccess::clearList()
{	n = 0;	}


#endif
