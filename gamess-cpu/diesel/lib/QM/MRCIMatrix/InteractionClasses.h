//***********************************************************************
//
//	Name:			InteractionClasses.h
//
//	Description:	
//
//	Author:			Michael Hanrath
//
//	Version:		0.0
//
//	Date:			31.08.1998
//
//
//
//
//
//***********************************************************************


#ifndef __InteractionClasses_h
#define __InteractionClasses_h

#include "../../../config.h"


template <class MatrixType, class VectorType> class HMatElements;
class BinomialCoefficient;
template <class TMOType> class	TableCase;
template <class TMOType> class	Configuration;

#include "../Configuration/DiffConf.h"
#include "../IntegralIndex/TwoElectronIntegralIndexUpdateMask.h"

template <class MatrixElementCalculator>
class InteractionClasses {
public:
	InteractionClasses(BinomialCoefficient *);
	~InteractionClasses();


	void	clear();

	

struct TCase {
	MatrixElementCalculator				*repMats;
	TableCase<MOType>					*tablecase;
	DiffConf<MOType>					diffConf;
	Configuration<MOType>				*to;
	TwoElectronIntegralIndexUpdateMask	updateMask;
	INT									extPosBOpen[2];
	INT									extPosBClosed[2];
	
	void	clearExtPosB();
	};
	
	
	TCase &	operator [] (INT i);

private:
static const INT NrCases = 16;
TCase cases[NrCases];

};



template <class MatrixElementCalculator>
inline
typename InteractionClasses<MatrixElementCalculator>::TCase &	InteractionClasses<MatrixElementCalculator>::operator [] (INT i)
{
	return	cases[i];
}


template <class MatrixElementCalculator>
inline
void	InteractionClasses<MatrixElementCalculator>::TCase::clearExtPosB()
{
	for ( INT i=0 ; i<2 ; i++ )
	{
		extPosBOpen[i] = -1;
		extPosBClosed[i] = -1;
	}
}

#endif

