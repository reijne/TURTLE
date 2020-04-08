//***********************************************************************
//
//	Name:			NTupelContainer.cc
//
//	Description:	base class for classes that contain
//					monades-, duades-, triade-integrals
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

#include "NTupelContainer.h"

#include "../Parallel/SharedMemory.h"

using namespace std;

NTupelContainer::NTupelContainer(
		MRMOs &mrmos,
		INT i1,
		const IrRep *irreps,
		INT i2)
{
	iIdx = NULL;
	p = NULL;
	n = NTupel = 0;
	
INT	nMO[4];

	map[0] = 0;
	map[1] = 1;
	map[2] = 2;

	for ( INT i=0 ; i<4 ; i++ )
		nMO[i] = mrmos.getNumberOfIrRepIntExt(irreps[i], i2 & (1<<i));

//	printf("Last:%d %d %d %d\n", nMO[3], nMO[2], nMO[1], nMO[0]);

	if ( nMO[0]*nMO[1]*nMO[2]*nMO[3]==0 )
		return;


	sub[0] = mrmos.getIrRepIntExtStart(irreps[3], i2 & 8);
	sub[1] = mrmos.getIrRepIntExtStart(irreps[2], i2 & 4);
	sub[2] = mrmos.getIrRepIntExtStart(irreps[1], i2 & 2);
	sub[3] = mrmos.getIrRepIntExtStart(irreps[0], i2 & 1);
//	printf("Sub:%d %d %d %d\n", sub[0], sub[1], sub[2], sub[3]);


	for ( INT i=0 ; i<3 ; i++ )
		if ( irreps[i]==irreps[i+1] && 
			(i2 & (1<<i)) > (i2 & (1<<(i+1))) )
			return;


INT	same[3];

	same[0] = (irreps[0]==irreps[1]) && (!(i2 & 1) == !(i2 & 2));
	same[1] = (irreps[1]==irreps[2]) && (!(i2 & 2) == !(i2 & 4));
	same[2] = (irreps[2]==irreps[3]) && (!(i2 & 4) == !(i2 & 8));


INT	sd[4];

	sd[3] = 0;
	for ( INT i=2 ; i>=0 ; i-- )
	{	if ( same[i] )
			sd[i] = sd[i+1] + (!(i1 & (1<<i)) ? 1 : 0);
		else
			sd[i] = 0;
		if ( nMO[i]-sd[i]<=0 )
			return;
	}



INT	*Idx[4];

	Idx[0] = new INT[nMO[0]];
	Idx[1] = new INT[nMO[1]];
	Idx[2] = new INT[nMO[2]];
	Idx[3] = new INT[nMO[3]];
INT m1=1;	
		fillIdx(nMO[0], Idx[0], m1);
INT m0=0;
	
INT	mult = nMO[0];
INT	k = 0;
	for ( INT i=0 ; i<3 ; i++ )
	{	if ( same[i] )
		{	if ( i1 & (1<<i) )
			{	fillIdx(nMO[i+1], Idx[i+1], m0);
//				printf("i1, i = %d %d\n", i1, i);
			}
			else
			{	fillIdx(nMO[i+1], Idx[i+1], Idx[k]);
				k = i+1;
//				printf("A i1, i = %d %d %d\n", i1, i, k);
			}
		}
		else
		{	mult = 1;
//			printf("mult1=%d\n", mult);
			for ( INT j=0 ; j<=i ; j++ )
				mult += Idx[j][nMO[j]-1-sd[j]];
		
//			printf("mult2=%d\n", mult);
			fillIdx(nMO[i+1], Idx[i+1], mult);

			k = i+1;
//			printf("B i1, i = %d %d %d\n", i1, i, k);
		}
	}

	delete[] Idx[0];
	iIdx = Idx[3];
	jIdx = Idx[2];
	kIdx = Idx[1];
//	printf("%d %d %d\n", iIdx[0], jIdx[0], kIdx[0]);
//	printf("%d %d %d\n", iIdx[nMO[3]-1], jIdx[nMO[2]-1], kIdx[nMO[1]-1]);
//	printf("%d %d %d\n", iIdx[nMO[3]-1-sd[3]], jIdx[nMO[2]-1-sd[2]], kIdx[nMO[1]-1-sd[1]]);
	n = iIdx[nMO[3]-1-sd[3]] + 
		jIdx[nMO[2]-1-sd[2]] + 
		kIdx[nMO[1]-1-sd[1]] + 
		nMO[0]-sd[0];

INT	permut[8] = {3, 2, 2, 1, 2, 2, 1, 1};
	NTupel = permut[i1];
	
//	cout << "n= " << n << endl;
//	cout << "nint= " << n*NTupel << endl;


	switch ( i1 ) {
	case 2:
		map[2] = 0;
		break;
	
	case 1:
	case 4:
	case 5:
		map[2] = 1;
		break;	
		
	case 3:
	case 6:
	case 7:
		map[1] = map[2] = 0;
		
	}
}


NTupelContainer::~NTupelContainer()
{
	if ( p )
		sharedMem->deallocate((void *) p);
	
	if ( iIdx )
	{	delete[] iIdx;
		delete[] jIdx;
		delete[] kIdx;
	}
}

void	NTupelContainer::setSharedMem(SharedMemory * _sharedMem)
{
	sharedMem = _sharedMem;
}


void	NTupelContainer::check()
{
	for ( INT i=0 ; i<n ; i++ )
		if ( p[3*i]!=777 || p[3*i+1]!=777 || p[3*i+2]!=777 )		
				cout << i << "   " << p[3*i] << "   " <<
				p[3*i+1] << "   " << p[3*i+2] << endl;
}

INT	NTupelContainer::allocate()
{
//	p = new IntegralType[n*NTupel];


//	p = new IntegralType[n*3];
	p = (IntegralType *) sharedMem->allocate(n*3, sizeof(IntegralType));
	memset(p, 0, n*3*sizeof(IntegralType));

//	printf("%x %d\n", p, n*NTupel);
//	p[0] = 7.44;
	if ( p )
		return 1;
	return 0;
}



void	NTupelContainer::fillIdx(INT n, INT *p, INT m)
{
	for ( INT i=0 ; i<n ; i++ )
		p[i] = i*m;
}


void	NTupelContainer::fillIdx(INT n, INT *p, INT *m)
{
	p[0] = 0;
	for ( INT i=1 ; i<n ; i++ )
		p[i] = p[i-1] + m[i-1];
}



