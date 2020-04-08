//***********************************************************************
//
//	Name:			MOAccess.cc
//
//	Description:	address lists for fast access
//
//	Author:			Michael Hanrath
//
//	Version:		0.0
//
//	Date:			26.09.1998
//
//
//
//
//
//***********************************************************************


#include "MOAccess.h"

#include "../MRMOs.h"
#include "MOIterator.h"


#include "../../MRTree/Diag/NExternalsDiag.h"


#include <string>



MOAccess::MOAccess(INT nTupel, const MRMOs *mrmos, INT maxInternalOpen,
	const NExternalsDiag &Psi0) :
	nTupel(nTupel), maxInternalOpen(maxInternalOpen), Psi0(Psi0)
{
	maxMO = mrmos->getMaxMO();
	
INT	n = 1;
	for ( INT i=0 ; i<nTupel ; i++ )
		n *= maxMO;
		
	pConf = new INT[n];
	
	pSAF = new INT * [maxInternalOpen+1];
	for ( INT i=0 ; i<=maxInternalOpen ; i++ )
	{
		pSAF[i] = new INT[n];
		memset(pSAF[i], 0, n*sizeof(INT));
	}

        INT* indSAF = new INT[maxInternalOpen+1];
	for ( IrRep irrep=0 ; irrep<mrmos->getNumberOfIrReps() ; irrep++ )
	{
	MOIterator	moIter(nTupel, mrmos, 1, irrep);
	INT	indConf = 0;

		for ( INT i=0 ; i<=maxInternalOpen ; i++ )
			indSAF[i] = 0;

		while ( !moIter.isEnd() )
		{
//			cout << "*** irrep=" << irrep << endl;
//			cout << moIter << endl;
		INT	m = 1;
		INT	k = 0;
		INT	ii = 1;
			for ( INT j=0 ; j<nTupel ; j++ )
			{
				k += m*(moIter.getMO(j)-1);
				m *= maxMO;
				if( j>0 )
				{
					if ( moIter.getMO(j)==moIter.getMO(j-1) )
						ii--;
					else
						ii++;
				}
//				cout << "k=" << k << ", " << moIter.getMO(j) << ", ii=" << ii << 
//				", ind[2]=" << ind[2] << ", " << Psi0.getNumberOfSpinAdaptedFunctions(
//								2 + ii) << endl;
			}
			
			for ( INT i=0 ; i<=maxInternalOpen ; i++ )
			{
				pSAF[i][k] = indSAF[i];
				indSAF[i] += Psi0.getNumberOfSpinAdaptedFunctions(
								i + ii);
			}

			pConf[k] = indConf++;
			

			moIter.next();
		}
	
			
	}
	delete[] indSAF;
	/*
	for ( INT i=0 ; i<=maxInternalOpen ; i++ )
	{
		cout << "open=" << i << endl;
		for ( INT j=1 ; j<=maxMO ; j++ )
		{
			for ( INT k=1 ; k<=maxMO ; k++ )
				cout << (*this)(i, j, k) << "\t";
			cout << endl;
		}
		cout << endl;
		cout << "-----------------------------" << endl;
		cout << endl;
	}
*/
}


MOAccess::~MOAccess()
{
	if ( pConf )
		delete pConf;

	if ( pSAF )
	{
		for ( INT i=0 ; i<=maxInternalOpen ; i++ )
			if ( pSAF[i] )
				delete pSAF[i];

		delete pSAF;
	}
}



