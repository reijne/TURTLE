//***********************************************************************
//
//	Name:			ExtPosContainer.cc
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




#include "ExtPosContainer.h"


MRMOs	ExtPosContainer::mrmos;

static const INT N = 16;

ExtPosContainer::ExtPosContainer(INT i1, const IrRep *irreps) :
	BaseContainer<NTupelContainer>(N)
{

	//	coding of i2:
	//
	//	i>=j>=k>=l
	//	i <--> internal  /  I <--> external
	//
	//			  bit 3	  bit 2	  bit 1	  bit 0
	//	(ij|kl)	 INT(i)  INT(j)  INT(k)  INT(l)		# internal
	//---------------------------------------------------------
	//	(ij|kl)		0		0		0		0			0

	//	(Ij|kl)		1		0		0		0			1
	//	(iJ|kl)		0		1		0		0			1
	//	(ij|Kl)		0		0		1		0			1
	//	(ij|kL)		0		0		0		1			1

	//	(IJ|kl)		1		1		0		0			2
	//	(Ij|Kl)		1		0		1		0			2
	//	(Ij|kL)		1		0		0		1			2
	//	(iJ|Kl)		0		1		1		0			2
	//	(iJ|kL)		0		1		0		1			2
	//	(ij|KL)		0		0		1		1			2

	//	(IJ|Kl)		1		1		1		0			3
	//	(IJ|kL)		1		1		0		1			3
	//	(Ij|KL)		1		0		1		1			3
	//	(iJ|KL)		0		1		1		1			3

	//	(IJ|KL)		1		1		1		1			4


	for ( INT i2=0 ; i2<N ; i2++ )
		if ( (!(i1 & 4) || !(i2 & 8) == !(i2 & 4)) &&
			 (!(i1 & 2) || !(i2 & 4) == !(i2 & 2)) &&
			 (!(i1 & 1) || !(i2 & 2) == !(i2 & 1)) )
		{//	printf("%d %d\n", i1, i2);
			p[i2] = new NTupelContainer(mrmos, i1, irreps, i2);
			if ( !p[i2]->getNumberOfIntegrals() )
			{	delete p[i2];
				p[i2] = NULL;
			}
		}
}


ExtPosContainer::~ExtPosContainer()
{
}


