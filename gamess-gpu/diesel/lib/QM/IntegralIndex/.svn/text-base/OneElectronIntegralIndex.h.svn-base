//***********************************************************************
//
//	Name:			OneElectronIntegralIndex.h
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

#ifndef __ONEELECTRONINTEGRALINDEX_H
#define __ONEELECTRONINTEGRALINDEX_H

#include "../../../config.h"

#include "../MO/MOType.h"

#include <string>
#include <iostream>
using std::ostream;

class OneElectronIntegralIndex {
public:
	OneElectronIntegralIndex();
	OneElectronIntegralIndex(INT i, INT j);
//	~OneElectronIntegralIndex();
	

	MOType *	getP();

	MOType &	operator [] (INT i);
	MOType	getIndex(INT i) const;

	MOType	getI() const;
	MOType	getJ() const;

	void	set(MOType i, MOType j);
	void	setI(MOType i);
	void	setJ(MOType i);



	friend ostream& operator<<(ostream& s, const OneElectronIntegralIndex &);

protected:
MOType	ind[2];
};



inline
OneElectronIntegralIndex::OneElectronIntegralIndex()
{
//	memset(ind, 0, 4*sizeof(INT));	
}


inline
OneElectronIntegralIndex::OneElectronIntegralIndex(
	INT _i, INT _j)
{
	ind[0] = _i;
	ind[1] = _j;
}

inline
MOType & 	OneElectronIntegralIndex::operator [] (INT i)
{	return ind[i];	}


inline
MOType *	OneElectronIntegralIndex::getP()
{	return	ind;	}

inline
MOType OneElectronIntegralIndex::getIndex(INT i) const
{	return ind[i];	}

inline
MOType	OneElectronIntegralIndex::getI() const
{	return ind[0];	}

inline
MOType	OneElectronIntegralIndex::getJ() const
{	return ind[1];	}


inline
void	OneElectronIntegralIndex::set(MOType i, MOType j)
{
	ind[0] = i;
	ind[1] = j;
}

inline
void	OneElectronIntegralIndex::setI(MOType i)
{	ind[0] = i;	}

inline
void	OneElectronIntegralIndex::setJ(MOType i)
{	ind[1] = i;	}




#endif
