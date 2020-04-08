//***********************************************************************
//
//	Name:			PropInts.cc
//
//	Description:	
//
//	Author:			Michael Hanrath
//
//	Version:		0.0
//
//	Date:			14.05.1998
//
//
//
//
//
//***********************************************************************


#include "PropInts.h"

#include "MOTrafo.h"
#include "DensityMatrix.h"

#include "../../Math/FortranNumeric/FortranLinkage.h"

#include <iostream>
#include <iomanip>

#include <string>

using namespace std;

extern "C" {
        void    FORTRAN_LINKAGE(propintr1)(
                INT *nReps
                );

        void    FORTRAN_LINKAGE(propintr2)(
                INT *nReps,
                INT *iReps
                );

        void    FORTRAN_LINKAGE(propintr3)(
                INT *nInts,
                INT *component,
   		INT *SymLbl
                );

        void    FORTRAN_LINKAGE(propintr4)(
                INT *mult,
                INT *nReps
                );

        void    FORTRAN_LINKAGE(propintr5)(
                INT *comp,
                INT *nInts,
                void *data
                );
};



char	*PropInts::OperatorNames[] = { "Mltpl  1", "Kinectic", "OneHam", "AngMom" };


PropInts::PropInts(Operator op, INT component, const char *filename)
{
	nInts = 0;
	p = NULL;

	{
	INT	iRc = -1;
	INT	iOpt = 0;
	INT	Lu = 12;	
//		FORTRAN_LINKAGE(opnone)(&iRc, &iOpt, filename, &Lu, strlen(filename));
//		cout << "iRc=" << iRc << endl;
	}	

	{
	INT	iRc = -1;
	INT	iOpt = 0;
	INT	iComp = 0;
	INT	iSymLbl = 1;
//		FORTRAN_LINKAGE(rdone)(&iRc, &iOpt, "nSym",
//			&iComp, &IrReps, &iSymLbl, strlen("nSym"));
		FORTRAN_LINKAGE(propintr1)(&IrReps);
//		cout << "IrReps=" << IrReps << endl;
	}	
	
	{
	INT	iRc = -1;
	INT	iOpt = 0;
	INT	iComp = 0;
	INT	iSymLbl = 1;

		inIrRep = new INT[IrReps];

		FORTRAN_LINKAGE(propintr2)(&IrReps,inIrRep);
//		FORTRAN_LINKAGE(rdone)(&iRc, &iOpt, "nBas",
//			&iComp, inIrRep, &iSymLbl, strlen("nBas"));
//		cout << inIrRep[0] << " " << inIrRep[1] << " "  << inIrRep[2] << " "  << inIrRep[3] << " " << endl;
	}	
	
	{
	INT	iRc = -1;
	INT	iOpt = 1;
	INT	iSymLbl = 1;
		FORTRAN_LINKAGE(propintr3)(&nInts,&component,&SymLbl);
//		FORTRAN_LINKAGE(rdone)(&iRc, &iOpt, OperatorNames[op],
//			&component, &nInts, &iSymLbl, strlen(OperatorNames[op]));
//		cout << "nInts=" << nInts << endl;
	}	
	
//	SymLbl = 1;
	{
	INT	iRc = -1;
	INT	iOpt = 6;

		p = new Type[nInts];

		FORTRAN_LINKAGE(propintr5)(&component,&nInts,p);
//		FORTRAN_LINKAGE(rdone)(&iRc, &iOpt, OperatorNames[op],
//			&component, p, &SymLbl, strlen(OperatorNames[op]));
	}	
	
	{
	INT	iRc = -1;
	INT	iOpt = 0;

//		FORTRAN_LINKAGE(clsone)(&iRc, &iOpt);
	}	

	init();


INT	i=0;	
	for ( INT iSym=0 ; iSym<IrReps ; iSym++ )
		for ( INT jSym=0 ; jSym<=iSym ; jSym++ )
			if ( (1 << (iSym ^ jSym)) & SymLbl )
				for ( INT iBas=0 ; iBas<inIrRep[jSym] ; iBas++ )
					for ( INT jBas=0 ; jBas<((iSym==jSym) ? iBas+1 : inIrRep[iSym]) ; jBas++ )
					{
//						cout << iSym << " " << jSym << " " << jBas << " " << iBas << " " << p[i] << endl;
						i++;
					}
						
/*
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	BUG in Beispielprogrammfragment von Fuelscher!!!!
	Schleifenindizes vertauscht!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!	
	for ( INT iSym=0 ; iSym<IrReps ; iSym++ )
		for ( INT jSym=0 ; jSym<=iSym ; jSym++ )
			if ( (1 << (iSym ^ jSym)) & SymLbl )
				for ( INT iBas=0 ; iBas<inIrRep[iSym] ; iBas++ )
					for ( INT jBas=0 ; jBas<((iSym==jSym) ? iBas+1 : inIrRep[jSym]) ; jBas++ )
						cout << iSym << " " << jSym << " " << iBas << " " << jBas << " " << p[i++] << endl;
*/						

	if ( i!=nInts )
	{
		cerr << "i=" << i << ", nInts=" << nInts << endl;
		cerr << "wrong number of integrals read. aborting" << endl;
		exit(1);
	}
}



PropInts::~PropInts()
{
	if ( p )
		delete p;
}



void	PropInts::transform(const MOTrafo &moTrafo)
{
/*	transformation of property integrals:
//
//         -----
//	       \
// 	|l>  =  >    a  |j>
//         /      j
//         -----
// 
// 
//  
//               -----                -----      
//               \                    \             
// 	<k|Op|l> = <  >    a   <i| | Op |  >    a   |j>  >
//               /      ki            /      lj      
//               -----                -----           
// 
//               -----      -----      
//               \          \             
// 	         = <  >    a     >    a     <i|Op|j> 
//               /      ki  /      lj      
//               -----      -----           
// 
// 	
// 	
*/

Type	*t1 = expand();
Type	*t2 = new Type[maxMO*maxMO];


	//	transform bra
	memset(t2, 0, maxMO*maxMO*sizeof(Type));
	for ( INT i=0 ; i<maxMO ; i++ )
		for ( INT j=0 ; j<maxMO ; j++ )
			for ( INT k=0 ; k<maxMO ; k++ )
				t2[i*maxMO + j] += t1[i*maxMO + k] * moTrafo(j+1, k+1);

	
/*
	collect(t2);
	cout << "________________________________________________" << endl;
	cout << *this << endl;
*/

	//	transform ket
	memset(t1, 0, maxMO*maxMO*sizeof(Type));
	for ( INT i=0 ; i<maxMO ; i++ )
		for ( INT j=0 ; j<maxMO ; j++ )
			for ( INT k=0 ; k<maxMO ; k++ )
				t1[i*maxMO + j] += t2[k*maxMO + i] * moTrafo(j+1, k+1);

	
	
	collect(t1);
	
	delete [] t1;
	delete [] t2;
}



PropInts::Type	PropInts::multDensity(const DensityMatrix &densityMatrix, bool anti) const
{
Type	v = 0;
INT	i=0;	
	for ( INT iSym=0 ; iSym<IrReps ; iSym++ )
		for ( INT jSym=0 ; jSym<=iSym ; jSym++ )
			if ( (1 << (iSym ^ jSym)) & SymLbl )
				for ( INT iBas=0 ; iBas<inIrRep[jSym] ; iBas++ )
					for ( INT jBas=0 ; jBas<((iSym==jSym) ? iBas+1 : inIrRep[iSym]) ; jBas++ )
					{
//						cout << iSym << " " << jSym << " " << jBas << " " << iBas << endl;
//						cout << p[i] << " " << densityMatrix(jSym, iSym, iBas+1, jBas+1) << " " << v << endl;
//							<< " " << densityMatrix(iSym, jSym, jBas+1, iBas+1) << endl;
						v += p[i]*densityMatrix(jSym, iSym, iBas+1, jBas+1);
						if ( fabs(densityMatrix(jSym, iSym, iBas+1, jBas+1))>0 )
						{
//							cout << jSym << " " << iSym << " " << iBas << " " << jBas << "    ";
//							cout << p[i] << " " << densityMatrix(jSym, iSym, iBas+1, jBas+1) << " " << v << endl;
						}
						if ( iSym!=jSym || iBas!=jBas )
						{
							v += (anti ? -1 : 1) * p[i]*densityMatrix(iSym, jSym, jBas+1, iBas+1);
							if ( fabs(densityMatrix(iSym, jSym, jBas+1, iBas+1))>0 )
							{
//								cout << iSym << " " << jSym << " " << jBas << " " << iBas << "    ";
//								cout << p[i] << " " << densityMatrix(iSym, jSym, jBas+1, iBas+1) << " " << v << endl;
							}
						}
						i++;
//						cout << v << endl << endl;
					}

	return v;
}

ostream & operator << (ostream &os, const PropInts &p)
{
PropInts::Type	*pp = p.expand();

	for ( INT i=0 ; i<p.maxMO ; i++ )
	{
		for ( INT j=0 ; j<p.maxMO ; j++ )
			os << setw(12) << pp[i*p.maxMO + j] << " ";
		os << endl;
	}
	
	delete pp;
	
	return os;
}

PropInts::Type *	PropInts::expand() const
{
Type	*pp = new Type[maxMO*maxMO];
	memset(pp, 0, maxMO*maxMO*sizeof(Type));
INT	i=0;	

	for ( INT iSym=0 ; iSym<IrReps ; iSym++ )
		for ( INT jSym=0 ; jSym<=iSym ; jSym++ )
			if ( (1 << (iSym ^ jSym)) & SymLbl )
				for ( INT iBas=0 ; iBas<inIrRep[jSym] ; iBas++ )
					for ( INT jBas=0 ; jBas<((iSym==jSym) ? iBas+1 : inIrRep[iSym]) ; jBas++ )
						pp[
							(getStartMO(iSym)+jBas-1)*maxMO +
							(getStartMO(jSym)+iBas-1)
							] =
						pp[
							(getStartMO(jSym)+iBas-1)*maxMO +
							(getStartMO(iSym)+jBas-1)
							] = p[i++];
	return pp;
}


void	PropInts::collect(const Type *pp) const
{
INT	i=0;	
	for ( INT iSym=0 ; iSym<IrReps ; iSym++ )
		for ( INT jSym=0 ; jSym<=iSym ; jSym++ )
			if ( (1 << (iSym ^ jSym)) & SymLbl )
				for ( INT iBas=0 ; iBas<inIrRep[jSym] ; iBas++ )
					for ( INT jBas=0 ; jBas<((iSym==jSym) ? iBas+1 : inIrRep[iSym]) ; jBas++ )
						p[i++] = pp[
							(getStartMO(iSym)+jBas-1)*maxMO +
							(getStartMO(jSym)+iBas-1)
							];
}
