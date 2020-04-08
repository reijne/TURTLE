//***********************************************************************
//
//	Name:			TwoElectronIntegralIndex.cc
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




#include "TwoElectronIntegralIndex.h"
#include "IndexMask.h"


template <class TMOType>
ostream& operator<<(ostream& s, const TwoElectronIntegralIndex<TMOType> & t)
{
	s << "(" << t.getI() << " " << t.getJ() << " | " <<
		t.getK() << " " << t.getL() << ")";
	return s;
}


/*
template <class TMOType>
TwoElectronIntegralIndex<TMOType>::TwoElectronIntegralIndex(
	const TwoElectronIntegralTriadeIndex & t)
{
}
*/

template <class TMOType>
void	TwoElectronIntegralIndex<TMOType>::setFromMask(IndexMask mask, MOType mo)
{
	for ( INT i=0 ; i<mask.getNumberOfFlags() ; i++ )
		if ( mask.getMask() & (1 << i) )
			ind[i] = mo;
}


template class TwoElectronIntegralIndex<MOType>;
template class TwoElectronIntegralIndex<GeneralizedMO>;

#define INSTANTIATE_TEMPLATE_FUNCTIONS \
template ostream& operator << (ostream& s, const TwoElectronIntegralIndex<T> & v);


#define T MOType
INSTANTIATE_TEMPLATE_FUNCTIONS
#undef T

#define T GeneralizedMO
INSTANTIATE_TEMPLATE_FUNCTIONS
#undef T
