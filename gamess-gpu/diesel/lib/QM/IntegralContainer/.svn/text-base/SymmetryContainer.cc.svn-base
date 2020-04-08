//***********************************************************************
//
//	Name:			SymmetryContainer.cc
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




#include "SymmetryContainer.h"

MRMOs	SymmetryContainer::mrmos;


static const INT N = 666;

SymmetryContainer::SymmetryContainer(INT i1) :
	BaseContainer<ExtPosContainer>(N)
{
ExtPosContainer::mrmos = mrmos;


//	printf("\n\n\n i1 = %d \n", i1);


	//	coding of i1:
	//
	//	i>=j>=k>=l
	//			  bit 2	  bit 1	  bit 0
	//	(ij|kl)	   i=j	   j=k	   k=l		dif. index	max irreps	# permut.
	//------------------------------------------------------------------------
	//	(ii|ii)		1		1		1			1			1			1
	//
	//	(ij|jj)		0		1		1			2			2			1
	//	(ii|ij)		1		1		0			2			2			1
	//
	//	(ii|jj)		1		0		1			2			2			2
	//
	//	(ij|kk)		0		0		1			3			2			2
	//	(ij|jk)		0		1		0			3			2			2
	//	(ii|jk)		1		0		0			3			2			2
	//
	//	(ij|kl)		0		0		0			4			4			3


IrRep	irreps[4];

	for ( irreps[3]=0 ; 
			irreps[3]<mrmos.getNumberOfIrReps() ; 
			irreps[3]++ )
		for ( irreps[2]=((i1 & 4) ? irreps[3] : 0) ;
				irreps[2]<=irreps[3] ;
				irreps[2]++ )
			for ( irreps[1]=((i1 & 2) ? irreps[2] : 0) ;
					irreps[1]<=((i1 & 2) ? irreps[2] : irreps[3]) ;
					irreps[1]++ )
				for ( irreps[0]=((i1 & 1) ? irreps[1] : 0) ;
					irreps[0]<=((i1 & 1) ? irreps[1] : 
						((irreps[3]==irreps[1]) ? irreps[2] : irreps[1])) ;
					irreps[0]++ )
				{//	printf("(%d %d|%d %d)\n", irreps[3], irreps[2], irreps[1], irreps[0]);
					if ( mrmos.getProd(
							mrmos.getProd(irreps[0], irreps[1]),
							mrmos.getProd(irreps[2], irreps[3])) == 0 
							&& irreps[2]>=irreps[1]
							)
						p[iM1over2(iM1over2(irreps[3]+1)+irreps[2]+1)+iM1over2(irreps[1]+1)+irreps[0]] =
							new ExtPosContainer(i1, irreps);
				}
}


SymmetryContainer::~SymmetryContainer()
{
}


