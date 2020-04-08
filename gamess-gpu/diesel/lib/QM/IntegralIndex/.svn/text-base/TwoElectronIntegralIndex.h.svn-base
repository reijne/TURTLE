//***********************************************************************
//
//	Name:			TwoElectronIntegralIndex.h
//
//	Description:	
//
//	Author:			Michael Hanrath
//
//	Version:		0.0
//
//	Date:			12. Feb 1998
//
//***********************************************************************

#ifndef __TWOELECTRONINTEGRALINDEX_H
#define __TWOELECTRONINTEGRALINDEX_H

#include "../../../config.h"

#include "../MO/MOType.h"
#include "../MO/GeneralizedMO.h"

#include <string>
#include <iostream>


class TwoElectronIntegralTriadeIndex;
class IndexMask;

template <class TMOType> class TwoElectronIntegralIndex;
template <class TMOType> ostream& operator<< (ostream& s, const TwoElectronIntegralIndex<TMOType> &);

template <class TMOType>
class TwoElectronIntegralIndex {
public:
	TwoElectronIntegralIndex();
	TwoElectronIntegralIndex(INT i, INT j, INT k, INT l);
//	~TwoElectronIntegralIndex();
	
	TwoElectronIntegralIndex(const TwoElectronIntegralTriadeIndex &);

	TMOType *	getP();

	TMOType &	operator [] (INT i);
	TMOType	getIndex(INT i) const;

	TMOType	getI() const;
	TMOType	getJ() const;
	TMOType	getK() const;
	TMOType	getL() const;

	void	set(TMOType i, TMOType j, TMOType k, TMOType l);
	void	setI(TMOType i);
	void	setJ(TMOType i);
	void	setK(TMOType i);
	void	setL(TMOType i);


	void	setFromMask(IndexMask mask, MOType mo);

	friend ostream& operator<< <TMOType> (ostream& s, const TwoElectronIntegralIndex<TMOType> &);

protected:
TMOType	ind[4];
};

template <class TMOType>	ostream& operator<<(ostream& s, const TwoElectronIntegralIndex<TMOType> &);


template <class TMOType>
inline
TwoElectronIntegralIndex<TMOType>::TwoElectronIntegralIndex()
{
//	memset(ind, 0, 4*sizeof(INT));	
}


template <class TMOType>
inline
TwoElectronIntegralIndex<TMOType>::TwoElectronIntegralIndex(
	INT _i, INT _j, INT _k, INT _l)
{
	ind[0] = _i;
	ind[1] = _j;
	ind[2] = _k;
	ind[3] = _l;
}

template <class TMOType>
inline
TMOType & 	TwoElectronIntegralIndex<TMOType>::operator [] (INT i)
{	return ind[i];	}


template <class TMOType>
inline
TMOType *	TwoElectronIntegralIndex<TMOType>::getP()
{	return	ind;	}

template <class TMOType>
inline
TMOType TwoElectronIntegralIndex<TMOType>::getIndex(INT i) const
{	return ind[i];	}

template <class TMOType>
inline
TMOType	TwoElectronIntegralIndex<TMOType>::getI() const
{	return ind[0];	}

template <class TMOType>
inline
TMOType	TwoElectronIntegralIndex<TMOType>::getJ() const
{	return ind[1];	}

template <class TMOType>
inline
TMOType	TwoElectronIntegralIndex<TMOType>::getK() const
{	return ind[2];	}

template <class TMOType>
inline
TMOType	TwoElectronIntegralIndex<TMOType>::getL() const
{	return ind[3];	}


template <class TMOType>
inline
void	TwoElectronIntegralIndex<TMOType>::set(
	TMOType i, TMOType j, TMOType k, TMOType l)
{
	ind[0] = i;
	ind[1] = j;
	ind[2] = k;
	ind[3] = l;
}

template <class TMOType>
inline
void	TwoElectronIntegralIndex<TMOType>::setI(TMOType i)
{	ind[0] = i;	}

template <class TMOType>
inline
void	TwoElectronIntegralIndex<TMOType>::setJ(TMOType i)
{	ind[1] = i;	}

template <class TMOType>
inline
void	TwoElectronIntegralIndex<TMOType>::setK(TMOType i)
{	ind[2] = i;	}

template <class TMOType>
inline
void	TwoElectronIntegralIndex<TMOType>::setL(TMOType i)
{	ind[3] = i;	}




#endif
