//***********************************************************************
//
//	Name:			TwoElectronIntegralIndexTemplate.h
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




#ifndef __TWOELECTRONINTEGRALINDEXTEMPLATE_H
#define __TWOELECTRONINTEGRALINDEXTEMPLATE_H

#include "../../../config.h"


#include "../MO/MOType.h"
#include "TwoElectronIntegralIndex.h"
#include "../MO/Iterators/MOListIterator.h"


class TwoElectronIntegralIndexTemplate :
	public TwoElectronIntegralIndex<MOType> {
public:
	TwoElectronIntegralIndexTemplate();
	TwoElectronIntegralIndexTemplate(INT i, INT j, INT k, INT l);

	// type conversion
	TwoElectronIntegralIndexTemplate(
		const TwoElectronIntegralIndex<GeneralizedMO> &);


	void	setFromTemplate(MOType *mo);
	void	setFromTemplate(MOType mo, INT n = 0);
	void	setFromTemplate(const MOListIterator &moiter);
	
	INT	getMaskI() const;
	INT	getMaskJ() const;
	INT	getMaskK() const;
	INT	getMaskL() const;

	void	setMaskI(INT flag);
	void	setMaskJ(INT flag);
	void	setMaskK(INT flag);
	void	setMaskL(INT flag);

	
private:
	void	calcStartEnd();
	
INT	mask[4];
INT	start;
INT	end;
};




inline
TwoElectronIntegralIndexTemplate::TwoElectronIntegralIndexTemplate()
{	
	memset(mask, -1, 4*sizeof(INT));
	start = 0;
	end = -1;
}

inline
TwoElectronIntegralIndexTemplate::TwoElectronIntegralIndexTemplate(
	INT i, INT j, INT k, INT l) :
	TwoElectronIntegralIndex<MOType>(i, j, k, l)
{	
	memset(mask, -1, 4*sizeof(INT));
	start = 0;
	end = -1;
}

inline
TwoElectronIntegralIndexTemplate::TwoElectronIntegralIndexTemplate(
	const TwoElectronIntegralIndex<GeneralizedMO> & twoeInd)
{
	start = 10;
	end = -10;
	for ( INT i=0 ; i<4 ; i++ )
		if ( twoeInd.getIndex(i).getType()==GeneralizedMO::Number )
		{	ind[i] = twoeInd.getIndex(i).getMONumber();
			mask[i] = -1;
		}
		else
		{	ind[i] = -1;
			mask[i] = twoeInd.getIndex(i).getSignature();
			if ( start>i )
				start = i;
			if ( i>end )
				end = i;
		}
}


inline
void	TwoElectronIntegralIndexTemplate::setFromTemplate(MOType *mo)
{
//	1. implementation
//	(highly efficient)
	for ( INT i=start ; i<=end ; i++ )
		if ( mask[i]>=0 )
			ind[i] = mo[mask[i]];

/*
//	2. implementation
//	(less efficient)
	for ( INT i=0 ; i<n ; i++ )
		ind[indind[i]] = mo[mask[i]];
*/
}



inline
void	TwoElectronIntegralIndexTemplate::setFromTemplate(MOType mo, INT n)
{
//	1. implementation
//	(highly efficient)
	for ( INT i=start ; i<=end ; i++ )
		if ( mask[i]==n )
			ind[i] = mo;
}

inline
void	TwoElectronIntegralIndexTemplate::setFromTemplate(
	const MOListIterator &moiter)
{
//	1. implementation
//	(highly efficient)
	for ( INT i=start ; i<=end ; i++ )
	{
		if ( mask[i]==0 )
			ind[i] = moiter.i;
		else
		if ( mask[i]==1 )
			ind[i] = moiter.j;
	}
}


inline
INT	TwoElectronIntegralIndexTemplate::getMaskI() const
{	return mask[0];	}

inline
INT	TwoElectronIntegralIndexTemplate::getMaskJ() const
{	return mask[1];	}

inline
INT	TwoElectronIntegralIndexTemplate::getMaskK() const
{	return mask[2];	}

inline
INT	TwoElectronIntegralIndexTemplate::getMaskL() const
{	return mask[3];	}


inline
void	TwoElectronIntegralIndexTemplate::calcStartEnd()
{
	start = 10;
	end = -10;
	for ( INT i=0 ; i<4 ; i++ )
	{	if ( start>i )
			start = i;
		if ( i>end )
			end = i;
	}
}

inline
void	TwoElectronIntegralIndexTemplate::setMaskI(INT flag)
{
	mask[0] = flag;
	calcStartEnd();
}

inline
void	TwoElectronIntegralIndexTemplate::setMaskJ(INT flag)
{
	mask[1] = flag;
	calcStartEnd();
}


inline
void	TwoElectronIntegralIndexTemplate::setMaskK(INT flag)
{
	mask[2] = flag;
	calcStartEnd();
}

inline
void	TwoElectronIntegralIndexTemplate::setMaskL(INT flag)
{
	mask[3] = flag;
	calcStartEnd();
}


#endif
