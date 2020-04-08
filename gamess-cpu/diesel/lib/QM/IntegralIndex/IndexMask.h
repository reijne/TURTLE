//***********************************************************************
//
//	Name:			IndexMask.h
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

#ifndef __INDEXMASK_H
#define __INDEXMASK_H

#include "../../../config.h"

#include "TwoElectronIntegralIndex.h"
#include "../MO/MOType.h"

#include <iostream>

typedef INT MaskType;

class IndexMask {
public:
	IndexMask(INT n);
	~IndexMask(){};


	void	setMask(TwoElectronIntegralIndex<MOType> ind, INT mo);

	void	setFlag(INT i, INT flag);

	INT	getFlag(INT i) const;
	MaskType	getMask() const;
	INT	getNumberOfFlags() const;


	friend ostream& operator<<(ostream& s, const IndexMask &);

private:
MaskType	mask;
INT	n;
};



inline
IndexMask::IndexMask(INT _n)
{	n = _n;	}


inline
void	IndexMask::setMask(TwoElectronIntegralIndex<MOType> ind, INT mo)
{
	mask = 0;
	for ( INT i=0 ; i<n ; i++ )
		mask |= (ind[i]==mo) << i;
}

inline
void	IndexMask::setFlag(INT i, INT flag)
{	mask = (mask & (~(1 << i))) | (flag << i);	}


inline
INT	IndexMask::getFlag(INT i) const
{	return (mask & (1 << i))!=0;	}

inline
MaskType	IndexMask::getMask() const
{	return mask;	}

inline
INT	IndexMask::getNumberOfFlags() const
{	return n;	}



#endif
