//***********************************************************************
//
//	Name:			TwoElectronIntegralTriadeIndex.h
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

#ifndef __TWOELECTRONINTEGRALTRIADEINDEX_H
#define __TWOELECTRONINTEGRALTRIADEINDEX_H

#include "../../../config.h"

#include "TwoElectronIntegralIndex.h"

#include <string>
#include <iostream>

class TwoElectronIntegralTriadeIndex : 
	public TwoElectronIntegralIndex<MOType> {
public:
	TwoElectronIntegralTriadeIndex();
	TwoElectronIntegralTriadeIndex(INT i, INT j, INT k, INT l, INT m);
//	~TwoElectronIntegralTriadeIndex();
	

	// type conversion
	TwoElectronIntegralTriadeIndex(const TwoElectronIntegralIndex<MOType> &);


	void	set(const TwoElectronIntegralIndex<MOType> &);
	
	INT	getM() const;
	
	void	setTriade();

	friend ostream& operator<<(ostream& s, const TwoElectronIntegralTriadeIndex &);



protected:
INT	m;
static unsigned short	IndexTable[64];
};


inline
TwoElectronIntegralTriadeIndex::TwoElectronIntegralTriadeIndex():
       TwoElectronIntegralIndex<MOType>(0,0,0,0)
{
//	m = 0;	
}


inline
TwoElectronIntegralTriadeIndex::TwoElectronIntegralTriadeIndex(
	INT _i, INT _j, INT _k, INT _l, INT _m) : 
	TwoElectronIntegralIndex<MOType>(_i, _j, _k, _l)
{	m = _m;	}


inline
INT	TwoElectronIntegralTriadeIndex::getM() const
{	return m;	}


inline
TwoElectronIntegralTriadeIndex::TwoElectronIntegralTriadeIndex(
	const TwoElectronIntegralIndex<MOType> & t)
{
INT	code = 	(t.getI()>t.getJ()) |
		(t.getI()>t.getK()) << 1 |
		(t.getI()>t.getL()) << 2 |
		(t.getJ()>t.getK()) << 3 |
		(t.getJ()>t.getL()) << 4 |
		(t.getK()>t.getL()) << 5;

INT	h = IndexTable[code];
	m = h & 3;
	code = h >> 8;
	
	ind[0] = t.getIndex((code >> 6) & 0x03);
	ind[1] = t.getIndex((code >> 4) & 0x03);
	ind[2] = t.getIndex((code >> 2) & 0x03);
	ind[3] = t.getIndex(code        & 0x03);
}


inline
void	TwoElectronIntegralTriadeIndex::set(
	const TwoElectronIntegralIndex<MOType> & t)
{
INT	code = 	(t.getI()>t.getJ()) |
		(t.getI()>t.getK()) << 1 |
		(t.getI()>t.getL()) << 2 |
		(t.getJ()>t.getK()) << 3 |
		(t.getJ()>t.getL()) << 4 |
		(t.getK()>t.getL()) << 5;

INT	h = IndexTable[code];
	m = h & 3;
	code = h >> 8;
	
	ind[0] = t.getIndex((code >> 6) & 0x03);
	ind[1] = t.getIndex((code >> 4) & 0x03);
	ind[2] = t.getIndex((code >> 2) & 0x03);
	ind[3] = t.getIndex(code        & 0x03);
}



inline
void	TwoElectronIntegralTriadeIndex::setTriade()
{
	if ( !m )
		return;

MOType	h;
	if ( m==1 )
	{	h = getL();
		setL(getK());
		setK(getJ());
		setJ(h);
	}
	else
	{	h = getK();
		setK(getJ());
		setJ(h);
	}
}

#endif
