//***********************************************************************
//
//	Name:			TwoElectronIntegralIndexUpdateMask.h
//
//	Description:	implements configuration handling
//					based on second quantization
//
//	Author:			Michael Hanrath
//
//	Version:		0.0
//
//	Date:			27.08.1998
//
//
//
//***********************************************************************

#ifndef __TwoElectronIntegralIndexUpdateMask_h
#define __TwoElectronIntegralIndexUpdateMask_h

#include "../../../config.h"



#include "TwoElectronIntegralIndex.h"
#include "../MO/MOType.h"



class TwoElectronIntegralIndexUpdateMask {
public:
	TwoElectronIntegralIndexUpdateMask();
	TwoElectronIntegralIndexUpdateMask(
		const TwoElectronIntegralIndex<MOType> &, const MOType *, INT n);
	~TwoElectronIntegralIndexUpdateMask() {};

	void	useOn(TwoElectronIntegralIndex<MOType> &, const MOType *) const;

	friend ostream & operator << (ostream &, const TwoElectronIntegralIndexUpdateMask &);
private:
INT	mask[4];
};


#include <string>


inline
TwoElectronIntegralIndexUpdateMask::TwoElectronIntegralIndexUpdateMask()
{
	for ( INT i=0 ; i<4 ; i++ )
		mask[i] = -1;
}



inline
TwoElectronIntegralIndexUpdateMask::TwoElectronIntegralIndexUpdateMask(
		const TwoElectronIntegralIndex<MOType> &index, const MOType *mo, INT n)
{
	for ( INT i=0 ; i<4 ; i++ )
	{
		mask[i] = -1;
		for ( INT j=0 ; j<n ; j++ )
			if ( ((TwoElectronIntegralIndex<MOType> &) index)[i]==mo[j] )
			{
				mask[i] = j;
				break;
			}
	}
}

inline
void	TwoElectronIntegralIndexUpdateMask::useOn(
	TwoElectronIntegralIndex<MOType> &index, const MOType *mo) const
{
	for ( INT i=0 ; i<4 ; i++ )
	{
		if ( mask[i]>=0 )
			index[i] = mo[mask[i]];
	}
}


#endif
