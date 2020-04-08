//***********************************************************************
//
//	Name:			TwoElectronIntegralCbExIndex.h
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

#ifndef __TWOELECTRONINTEGRALCBEXINDEX_H
#define __TWOELECTRONINTEGRALCBEXINDEX_H

#include "../../../config.h"

#include "TwoElectronIntegralTriadeIndex.h"

#include <string>
#include <iostream>

class TwoElectronIntegralCbExIndex : 
	public TwoElectronIntegralTriadeIndex {
public:
	TwoElectronIntegralCbExIndex() {}
	TwoElectronIntegralCbExIndex(INT i, INT j, INT k, INT l, INT cb, INT ex);
//	~TwoElectronIntegralCbExIndex();
	

	// type conversion
	TwoElectronIntegralCbExIndex(
		const TwoElectronIntegralIndex<MOType> &,
		const TwoElectronIntegralIndex<MOType> &
		);
	
	void	set(const TwoElectronIntegralIndex<MOType> &ind)
	{	TwoElectronIntegralTriadeIndex::set(ind);	}

	void	set(
		const TwoElectronIntegralIndex<MOType> &,
		const TwoElectronIntegralIndex<MOType> &
	);

	INT	getCoulombTriade() const;
	INT	getExchangeTriade() const;

	friend ostream& operator<<(ostream& s, const TwoElectronIntegralCbExIndex &);



private:
INT	n;
};



inline
TwoElectronIntegralCbExIndex::TwoElectronIntegralCbExIndex(
	INT _i, INT _j, INT _k, INT _l, INT _cb, INT _ex) : 
	TwoElectronIntegralTriadeIndex(_i, _j, _k, _l, _cb)
{
	n = _ex;
}


inline
INT	TwoElectronIntegralCbExIndex::getCoulombTriade() const
{	return m;	}

inline
INT	TwoElectronIntegralCbExIndex::getExchangeTriade() const
{	return n;	}


inline
TwoElectronIntegralCbExIndex::TwoElectronIntegralCbExIndex(
	const TwoElectronIntegralIndex<MOType> & cb,
	const TwoElectronIntegralIndex<MOType> & ex)
{
INT	code = 	(cb.getI()>cb.getJ()) |
		(cb.getI()>cb.getK()) << 1 |
		(cb.getI()>cb.getL()) << 2 |
		(cb.getJ()>cb.getK()) << 3 |
		(cb.getJ()>cb.getL()) << 4 |
		(cb.getK()>cb.getL()) << 5;

INT	h = IndexTable[code];
	m = h & 3;
	code = h >> 8;
	n = (h >> 2) & 0x3f;
	h = (cb.getL()==ex.getL()) | (cb.getK()==ex.getK()) << 1;
	n = (n >> (2*h)) & 3;

	ind[0] = cb.getIndex((code >> 6) & 0x03);
	ind[1] = cb.getIndex((code >> 4) & 0x03);
	ind[2] = cb.getIndex((code >> 2) & 0x03);
	ind[3] = cb.getIndex(code        & 0x03);
}


inline
void	TwoElectronIntegralCbExIndex::set(
	const TwoElectronIntegralIndex<MOType> & cb,
	const TwoElectronIntegralIndex<MOType> & ex)
{
INT	code = 	(cb.getI()>cb.getJ()) |
		(cb.getI()>cb.getK()) << 1 |
		(cb.getI()>cb.getL()) << 2 |
		(cb.getJ()>cb.getK()) << 3 |
		(cb.getJ()>cb.getL()) << 4 |
		(cb.getK()>cb.getL()) << 5;

INT	h = IndexTable[code];
	m = h & 3;
	code = h >> 8;
	n = (h >> 2) & 0x3f;
	h = (cb.getL()==ex.getL()) | (cb.getK()==ex.getK()) << 1;
	n = (n >> (2*h)) & 3;

	ind[0] = cb.getIndex((code >> 6) & 0x03);
	ind[1] = cb.getIndex((code >> 4) & 0x03);
	ind[2] = cb.getIndex((code >> 2) & 0x03);
	ind[3] = cb.getIndex(code        & 0x03);
}

#endif
