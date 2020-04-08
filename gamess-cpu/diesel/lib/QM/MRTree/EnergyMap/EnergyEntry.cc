//***********************************************************************
//
//	Name:			EnergyEntry.cc
//
//	Description:	
//
//	Author:			Michael Hanrath
//
//	Version:		0.0
//
//	Date:			19.04.1997
//
//
//
//
//***********************************************************************

#include "EnergyEntry.h"

#include <iostream>
#include <string>
#include <math.h>

using namespace std;

EnergyEntry::Status	EnergyEntry::select(
	EnergyType threshold,
	SelectionMode mode)
{
INT	Degen = 0;
EnergyEntry::Status	flag = NotSelected;

	dE = new EnergyType[nRoots];
	memset(dE, 0, nRoots*sizeof(EnergyType));

	ci = new PTCIType[nRoots*nSAFs];
PTCIType	*pci = ci;
	
//
//                            Main  SAF_Main(k)                          2
//                          | -----    -----                            |
//                          | \        \                                |
//                          |  >        >    c_klr <phi_kl | H | phi_ij>|
//           Gen  SAF_Gen(i)| /        /                                |
//          -----   -----   | -----    -----                            |  
//          \       \          k        l
// eps_ges = >       >    ------------------------------------------------
//          /       /
//          -----   -----           E_0r - <phi_ij | H | phi_ij>
//            i       j
//
//                   /|\ 
//                    |
//                    |    
//                    +----- this routine calculates this part of the sums
//
//	r: root number

	// configurations are selected if threshold is reached for any root
	switch ( mode ) {
	case AmountSquareOfSum:
		{
		EnergyType	*pnom = nominator;
		EnergyType	*pdenom = denominator;
//			cout << "~s ";
			for ( INT i=0 ; i<nRoots ; i++ )
			{
				for ( INT j=0 ; j<nSAFs ; j++ )
				{
				double	h;
//					cout << ":::" << *pnom << " " << *pdenom << endl;
					Degen |= (fabs(*pdenom) < DegenerationThreshold);
					h = *pnom / *pdenom++;
					*pci++ = h;
					dE[i] += *pnom++ * h;
//					dE[i] += h; pnom++;
				}
//				cout << "i=" << i << ", dE[i]=" << dE[i] << endl;
//				cout << " " << dE[i];

//	cout << dE[i] << " " << nSAFs << endl;
//				if ( fabs(dE[i])>=threshold )
				if ( fabs(dE[i]/nSAFs)>=threshold )
					flag = EnergyEntry::Selected; 
			}
		}
		break;

	case SumOfAmountSquares:
		{
			return NotSelected;
		}
		break;
	}
	delete[] nominator;
	nominator = NULL;
	delete[] denominator;
	denominator = NULL;
	if ( Degen )
		return SelectionFlagserated;
	else
		return flag;
}

const EnergyType EnergyEntry::DegenerationThreshold = 1e-2;
