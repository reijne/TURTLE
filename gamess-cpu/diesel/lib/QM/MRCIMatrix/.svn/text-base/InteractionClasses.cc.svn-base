//***********************************************************************
//
//	Name:			InteractionClasses.cc
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




#include "InteractionClasses.h"

#include "../Configuration/Configuration.h"
#include "../Configuration/TableCase.h"



template <class MatrixElementCalculator>
InteractionClasses<MatrixElementCalculator>::InteractionClasses(BinomialCoefficient *binom)
{
	for ( INT i=0 ; i<NrCases ; i++ )
	{
		cases[i].tablecase = new TableCase<MOType>(binom);
		cases[i].to = (Configuration<MOType> *) &cases[i].diffConf.getTo();
	}
	clear();
}


template <class MatrixElementCalculator>
InteractionClasses<MatrixElementCalculator>::~InteractionClasses()
{
	for ( INT i=0 ; i<NrCases ; i++ )
		delete cases[i].tablecase;
}


template <class MatrixElementCalculator>
void	InteractionClasses<MatrixElementCalculator>::clear()
{
	for ( INT i=0 ; i<NrCases ; i++ )
		cases[i].repMats = NULL;
}

#include "../RepresentationMatrices/HMatElements.h"
#include "../RepresentationMatrices/MRMPH0MatElements.h"

template class InteractionClasses<HMatElements<float, float> >;
template class InteractionClasses<HMatElements<double, float> >;
template class InteractionClasses<HMatElements<float, double> >;
template class InteractionClasses<HMatElements<double, double> >;

template class InteractionClasses<MRMPH0MatElements<float, float> >;
template class InteractionClasses<MRMPH0MatElements<double, float> >;
template class InteractionClasses<MRMPH0MatElements<float, double> >;
template class InteractionClasses<MRMPH0MatElements<double, double> >;
