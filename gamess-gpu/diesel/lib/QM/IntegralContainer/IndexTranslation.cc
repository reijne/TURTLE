//***********************************************************************
//
//	Name:			IndexTranslation.cc
//
//	Description:	
//
//	Author:			Michael Hanrath
//
//	Version:		0.0
//
//	Date:			03.09.1996
//
//
//
//
//
//***********************************************************************




#include "IndexTranslation.h"


IndexTranslation::IndexTranslation(INT maxMO)
{
	for ( INT i=0 ; i<5 ; i++ )
		Index[i] = new INT[maxMO+2];


	for ( INT i=0 ; i<maxMO+2 ; i++ )
	{	Index[0][i] = 0;
		Index[1][i] = i;
		Index[2][i] = i*(i-1)/2;
		Index[3][i] = i*(i-1)*(i-2)/6;
		Index[4][i] = (((i*(i-1)/2)*(i-2)/3)*(i-3)/4);
	}
}



IndexTranslation::~IndexTranslation()
{
	for ( INT i=0 ; i<5 ; i++ )
		delete[] Index[i];
}
