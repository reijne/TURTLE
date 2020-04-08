//***********************************************************************
//
//	Name:			RepDiag.cc
//
//	Description:	stores information for an efficient calculation
//					of P=5 diagonal representation matrices
//
//					
//
//	Author:			Michael Hanrath
//
//	Version:		0.0
//
//	Date:			08.04.1997
//
//
//
//
//
//***********************************************************************




#include "RepDiag.h"



template <class RepMatType, class VectorType>
RepDiag<RepMatType, VectorType>::RepDiag(
	const HMatElements<RepMatType, VectorType>	*const *_repMatsP5,
	Configuration<MOType> conf,
	INT minOpenShells,
	INT maxOpenShells)
{
	repMatsP5 = _repMatsP5;
	intConf = conf;
	DiagIntPart = new RepMatType[maxOpenShells+1];
	for ( INT i=minOpenShells ; i<=maxOpenShells ; i+=2 )
		DiagIntPart[i] = repMatsP5[i]->getCase14DiagTeil(conf);
		
	DiagExtPart = NULL;
	DiagExtPartLast = NULL;
	diag = NULL;
	pP5 = NULL;
	pP5Last = NULL;
}



template <class RepMatType, class VectorType>
RepDiag<RepMatType, VectorType>::~RepDiag()
{
	delete[] DiagIntPart;
	
	if ( DiagExtPart )
	{
		delete[] DiagExtPart;
		delete[] diag;
	}

	if ( pP5 )
		delete[] pP5;

	if ( pP5Last )
		delete[] pP5Last;
		
	if ( DiagExtPartLast )
		delete[] DiagExtPartLast;
}




template <class RepMatType, class VectorType>
void	RepDiag<RepMatType, VectorType>::addExt(
	Configuration<MOType> conf,
	INT totalOpenShells,
	INT	moExtPos1,
	INT moExtPos2)
{
//	cout << conf << endl;
//	printf("%x %x %x %x %x\n", DiagExtPart, diag, pP5, pP5Last, DiagExtPartLast);
	if ( DiagExtPart )
	{
		delete[] DiagExtPart;
		delete[] diag;
	}

	if ( pP5 )
		delete[] pP5;

	if ( pP5Last )
		delete[] pP5Last;
		
	if ( DiagExtPartLast )
		delete[] DiagExtPartLast;

//	cout << *repMatsP5[totalOpenShells] << endl;

	dim = repMatsP5[totalOpenShells]->getNumberOfRows();
	DiagExtPart = new RepMatType[dim];
	diag = new RepMatType[dim];
	restConf = intConf;
	restConf.append(conf);
//	cout << "intConf=" << intConf << endl;
//	cout << "restConf=" << restConf << endl;
	for ( INT i=0 ; i<dim ; i++ )
		DiagExtPart[i] = DiagIntPart[totalOpenShells];

	nOpen = totalOpenShells - restConf.getNumberOfOpenShells();
	repMatP5 = repMatsP5[totalOpenShells];

//	printf("A %d %d %d     %d\n", dim, totalOpenShells, nOpen, restConf.getNumberOfOpenShells()*dim);
	DiagExtPartLast = NULL;
	pP5 = NULL;
	pP5Last = NULL;

	switch ( nOpen ) {
	case 0:
		moExtPos2 = moExtPos1 = -1;
		break;
		
	case 1:
		if ( restConf.getNumberOfOpenShells() )
			pP5 = new RepMatType[restConf.getNumberOfOpenShells()*dim];
		moExtPos2 = -1;
		break;
	
	case 2:
		if ( moExtPos1>moExtPos2 )
		{
		INT	h = moExtPos1;
			moExtPos1 = moExtPos2;
			moExtPos2 = h;
		}
		if ( restConf.getNumberOfOpenShells() )
			pP5 = new RepMatType[restConf.getNumberOfOpenShells()*dim];
		pP5Last = new RepMatType[(restConf.getNumberOfOpenShells()+1)*dim];
		DiagExtPartLast = new RepMatType[dim];
		break;
	}	

INT	l = 0;
MOType	moi = 0, moj = 0;
TwoElectronIntegralIndex<MOType>	ind1;
TwoElectronIntegralTriadeIndex	triadeIndex;
IntegralType	intBuf[MAXOPENSHELLS*(MAXOPENSHELLS+1)/2];
INT	ext[MAXOPENSHELLS*(MAXOPENSHELLS+1)/2];

INT	isExti;
INT	subi = 0;

//	cout << restConf.getNumberOfOpenShells() << endl;
//	cout << "totalOpenShells=" << totalOpenShells << " nOpen=" << nOpen << " moExtPos1=" << moExtPos1 << " moExtPos2=" << moExtPos2 << endl;
	
	// prepare diagonals of P=5 matrices for an easy later processing
	// this is a bit tricky ...
	
	if ( 0==moExtPos1 || 0==moExtPos2 )
		subi++;
	for ( INT i=1 ; i<totalOpenShells ; i++ )
	{
		isExti = 0;
		if ( i==moExtPos1 )
		{
			isExti = 1;
			subi++;
		}
		else
		if ( i==moExtPos2 )
		{
			isExti = 2;
			subi++;
		}
		else
			moi = restConf.getOpenShell(i - subi);

	INT	subj = 0;
		for ( INT j=0 ; j<i ; j++ )
		{
		INT	isExtj = 0;
			if ( j==moExtPos1 )
			{
				isExtj = 1;
				subj++;
			}
			else
			if ( j==moExtPos2 )
			{
				isExtj = 2;
				subj++;
			}
//			cout << "A i=" << i << " j=" << j << " isExti" << isExti << " isExtj" << isExtj << endl;
//			cout << subi << " " << subj << endl;
			if ( !(isExtj || isExti) )
			{
				moj = restConf.getOpenShell(j - subj);
				ind1.set(moi, moj, moj, moi);
				triadeIndex.set(ind1);
//				cout << "^^^^^^^^^" << ind1 << endl;
//				intBuf[l] = (*HMatElements<RepMatType, VectorType>::iijj)[triadeIndex];
				intBuf[l] = (*HMatElements<RepMatType, VectorType>::int4Container)[triadeIndex];
			}
//			cout << "i=" << i << " j=" << j << " isExti" << isExti << " isExtj" << isExtj << endl;
			ext[l++] = isExti + isExtj;
		}
	}

	const RepMatType	*p = 
		repMatsP5[totalOpenShells]->getCbMat()->getP();
	for ( INT i=0 ; i<dim ; i++ )
	{
	INT lll1 = 0;
	INT lll2 = 0;
		for ( INT ll=0 ; ll<l ; ll++ )
			switch ( ext[ll] ) {
			case 0:
				DiagExtPart[i] += intBuf[ll] * *p++;
				break;
				
			case 1:
//				cout << "lll1=" << lll1 << " " << *p << endl;
				pP5[lll1++ * dim + i] = *p++;
				break;
				
			case 2:
			case 3:
//				cout << "lll2=" << lll2 << " " << *p << endl;
				pP5Last[lll2++ * dim + i] = *p++;	
				break;
			}	
		// skip off-diagonal elements
		p += l*(i+1);
	}
}


template class RepDiag<float, float>;
template class RepDiag<double, float>;
template class RepDiag<float, double>;
template class RepDiag<double, double>;
