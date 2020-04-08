//***********************************************************************
//
//	Name:			RepresentationMatrixFortranInterface.cc
//
//	Description:	interface to the implementation in Fortran 
//					by Volker Pleﬂ
//					does the neccesary initializations and
//					calculates the representation matrices of the 
//					symmetric group approach (SGA) for one case
//					(matrices: "umat", "umat1", "emat", "hp5dar")
//
//					
//
//	Author:			Michael Hanrath
//
//	Version:		0.0
//
//	Date:			27.08.1996
//
//
//
//
//
//***********************************************************************

#include "RepresentationMatrixFortranInterface.h"
#include "../Configuration/ConfigurationGlobals.h"
#include "RepresentationMatrices.h"

#include "../../Math/FortranNumeric/FortranLinkage.h"

#include <stdio.h>
#include <stdlib.h>
#include <string>

#include "../IO/Verbosity.h"



template <class RepMatType>
SpinEigenFunctionDegeneration *	
	RepresentationMatrixFortranInterface<RepMatType>::spinEigs = NULL;

extern "C" {
	void	FORTRAN_LINKAGE(buminit)(
		INT *multiplicity,
		INT *matDim,
		INT *maxOp,
		double **pUMat,
		double **pUMat1,
		double **pEMat,
//		double **php5dar,
		INT	*verbose
		);
		
	void	FORTRAN_LINKAGE(calcbum)(
		INT *idsk, INT *ipfa, INT *irfa1,
		INT *iqrfa,INT *iqlfa, 
		INT *kml, INT *kmj, 
		INT *inopen, INT *ibob
		);

	void	FORTRAN_LINKAGE(calcdarhp5d)(
		INT *idsk, INT *kml, INT *inopen, double* pp
		);

	void	FORTRAN_LINKAGE(calcdarhp5f)(
		INT *idsk, INT *kml, INT *inopen, float* pp
		);
}

template <class RepMatType>
RepresentationMatrixFortranInterface<RepMatType>::RepresentationMatrixFortranInterface(
	INT	multiplicity)
{
	spinEigs = new SpinEigenFunctionDegeneration(multiplicity, MAXOPENSHELLS);

	if ( verbosity.isActive(Verbosity::SGA) )
		cout << *spinEigs << endl;

INT	verbose = verbosity.isActive(Verbosity::SGA);
	FORTRAN_LINKAGE(buminit)(&multiplicity, &matDim, &maxOp,
		&pUMat, &pUMat1, &pEMat, &verbose);
//		&pUMat, &pUMat1, &pEMat, &php5dar, &verbose);
		
/*	printf("matDim = %d\n", matDim);
	printf("Umat[0] = %lf\n", pUMat[0]);
	printf("Umat[matDim] = %lf\n", pUMat[matDim]);
	printf("Umat[1] = %lf\n", pUMat[1]);
*/

	RepresentationMatrices<RepMatType>::repMatFInt = this;	
}


template <class RepMatType>
RepresentationMatrixFortranInterface<RepMatType>::~RepresentationMatrixFortranInterface()
{
	delete spinEigs;
}


template <class RepMatType>
void	RepresentationMatrixFortranInterface<RepMatType>::calc(
	RepMatType *p1, RepMatType *p2,
	INT dK, INT P, INT R, INT qR, INT qL,
	INT kml, INT kmj, INT numberOfOpenShells
	)
{
INT	absdK = abs(dK);
INT	ibob;
RepMatType	*p;
INT	flag, eflag;

//FD also done in calcdarhp5d
//FD	memset(pUMat, 0, matDim*matDim*sizeof(RepMatType));


//	printf("%d %d %d %d %d %d %d %d\n", dK, P, R, qR, qL,
//	kml, kmj, numberOfOpenShells);

	if ( P==5 )
	{	p = p1;
		INT	n = numberOfOpenShells*(numberOfOpenShells-1)/2;
      		switch (sizeof(RepMatType))
		{
		case 8: FORTRAN_LINKAGE(calcdarhp5d)(&absdK, &kml, &numberOfOpenShells,(double*)p); break;
		case 4: FORTRAN_LINKAGE(calcdarhp5f)(&absdK, &kml, &numberOfOpenShells,(float*)p); break;
                default: printf ("remco is gek\n"); exit(1); break;
                }
//		p = p1;
//		INT	n = numberOfOpenShells*(numberOfOpenShells-1)/2;
//		cout << "*** P=5 *** n="<<n << endl;
//		for ( INT i=0 ; i<kml ; i++ )
//			for ( INT j=0 ; j<=i ; j++ )
//				for ( INT k=0 ; k<n ; k++ )
//					*p++ = php5dar[k + j*maxOp + i*maxOp*matDim];
		return;
	}
	else
	{
		if ( absdK==2 && R>=4 )
			R = 7-R;

		if ( dK==0 )
		{	ibob = 0;	// ibob ohne Bedeutung!!!
			FORTRAN_LINKAGE(calcbum)(
				&absdK, &P, &R, &qL, &qR, &kml, &kmj, &numberOfOpenShells, &ibob);
			flag = ibob;	// ibob wird von calcbum auf qR>qL gesetzt		
			eflag = 0;		// fuer emat laeuft immer zuerst qR, dann qL
		}
		else
		if ( dK>0 )	
		{	ibob = 0;	// Eingabeparameter
					// ibob=0 ==> niedrigeres l‰uft im EMat-Feld zuerst 
					// kein Einfluss auf UMat
					// (UMat immer aufrecht stehend)
			FORTRAN_LINKAGE(calcbum)(
				&absdK, &P, &R, &qR, &qL, &kmj, &kml, &numberOfOpenShells, &ibob);
//???????????????????????????????????????????????????????????????????????
			flag = 1;
//???????????????????????????????????????????????????????????????????????
//			flag = 0;
//???????????????????????????????????????????????????????????????????????
			eflag = 0;
		}
		else	
		{
			ibob = 0;
			FORTRAN_LINKAGE(calcbum)(
				&absdK, &P, &R, &qL, &qR, &kml, &kmj, &numberOfOpenShells, &ibob);

//???????????????????????????????????????????????????????????????????????
//			flag = 1;
//???????????????????????????????????????????????????????????????????????
			flag = 0;
//???????????????????????????????????????????????????????????????????????
			eflag = 1;
		}

/*		cout << "flag, ibob= " << flag << ", " << ibob << endl;

		for ( INT i=0 ; i<5 ; i++ )
		{	for ( INT j=0 ; j<5 ; j++ )
				cout << pUMat[i*matDim + j] << " ";
			cout << endl;
		}
*/

		if ( flag )
		{
			//	copy coulomb representation matrices
			p = p1;
			for ( INT i=0 ; i<kml ; i++ )
				for ( INT j=0 ; j<kmj ; j++ )
				{
//					cout << i << " " << j << " " << i*matDim + j << endl;
					*p++ = pUMat[i*matDim + j];
				}

			if ( P==1 )
			{
			//	copy exchange representation matrices
				p = p2;
				for ( INT i=0 ; i<kml ; i++ )
					for ( INT j=0 ; j<kmj ; j++ )
						*p++ = pUMat1[i*matDim + j];
			}
		}
		else
		{
			//	copy coulomb representation matrices
			p = p1;
			for ( INT i=0 ; i<kml ; i++ )
				for ( INT j=0 ; j<kmj ; j++ )
				{
//					cout << i << " " << j << " " << j*matDim + i << endl;
					*p++ = pUMat[j*matDim + i];
				}
			
			
			if ( P==1 )
			{
			//	copy exchange representation matrices
				p = p2;
				for ( INT i=0 ; i<kml ; i++ )
					for ( INT j=0 ; j<kmj ; j++ )
						*p++ = pUMat1[j*matDim + i];
			}
		}

		if ( P==3 )
		{	if ( eflag )
			{
			//	copy P=3 representation matrices
				p = p2;
			INT	n = numberOfOpenShells-absdK-1;
/*				cout << "----------------" << endl;
				for ( INT i=0 ; i<kml*kmj*n ; i++ )
						cout << pEMat[i] << endl;
				cout << "----------------" << endl;
*/				for ( INT i=0 ; i<kml ; i++ )
					for ( INT j=0 ; j<kmj ; j++ )
						for ( INT k=0 ; k<n ; k++ )
						{//	cout << pEMat[i*kmj*n + j*n + k] << endl;
							*p++ = pEMat[i*kmj*n + j*n + k];
						}
			}
			else
			{
			//	copy P=3 representation matrices
				p = p2;
			INT	n = numberOfOpenShells-absdK-1;
/*				cout << "----------------" << endl;
				for ( INT i=0 ; i<kml*kmj*n ; i++ )
						cout << pEMat[i] << endl;
				cout << "----------------" << endl;
*/				for ( INT i=0 ; i<kml ; i++ )
					for ( INT j=0 ; j<kmj ; j++ )
						for ( INT k=0 ; k<n ; k++ )
						{//	cout << pEMat[j*kml*n + i*n + k] << endl;
							*p++ = pEMat[j*kml*n + i*n + k];
						}
			}
		}
	}
/*	cout << "-------------------------------------------------------" << endl;
	p = p1;
	for ( INT i=0 ; i<kml ; i++ )
	{	for ( INT j=0 ; j<kmj ; j++ )
			cout << *p++ << " ";
		cout << endl;
	}
	cout << "=======================================================" << endl;
*/
}


template class RepresentationMatrixFortranInterface<float>;
template class RepresentationMatrixFortranInterface<double>;
