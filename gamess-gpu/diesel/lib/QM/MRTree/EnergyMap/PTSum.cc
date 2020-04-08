//***********************************************************************
//
//	Name:			PTSum.cc
//
//	Description:	sums up energy estimations from
//					perturbation theory
//
//	Author:			Michael Hanrath
//
//	Version:		0.0
//
//	Date:			13.05.1997
//
//
//
//***********************************************************************

#include "PTSum.h"

#include <string>
#include <iomanip>

using  namespace std;

PTSum::PTSum(EnergyType SelectionThreshold, INT _divideByCSFs, INT _nRoots)
{
	nThresh = 10;
	nRoots = _nRoots;
	divideByCSFs = _divideByCSFs;
	thresh = new EnergyType[nThresh];
	sum = new EnergyType[nThresh*nRoots];
	countConf = new INT[nThresh];
	countCSF = new INT[nThresh];

	thresh[3] = SelectionThreshold;
	for ( INT i=4 ; i<nThresh ; i++ )
		thresh[i] = thresh[i-1] + thresh[3]/2;
	for ( INT i=2 ; i>=0 ; i-- )
		thresh[i] = thresh[i+1]/2;
		
	memset(sum, 0, nThresh*nRoots*sizeof(EnergyType));
	memset(countConf, 0, nThresh*sizeof(INT));
	memset(countCSF, 0, nThresh*sizeof(INT));
}

PTSum::PTSum(INT _nThresh, const EnergyType *e, INT _divideByCSFs, INT _nRoots)
{
	nThresh = _nThresh;
	nRoots = _nRoots;
	divideByCSFs = _divideByCSFs;
	thresh = new EnergyType[nThresh];
	sum = new EnergyType[nThresh*nRoots];
	countConf = new INT[nThresh];
	countCSF = new INT[nThresh];

	memcpy(thresh, e, nThresh*sizeof(EnergyType));	
	memset(sum, 0, nThresh*nRoots*sizeof(EnergyType));
	memset(countConf, 0, nThresh*sizeof(INT));
	memset(countCSF, 0, nThresh*sizeof(INT));
}

PTSum::PTSum(INT _nThresh, INT _divideByCSFs, INT _nRoots)
{
	nThresh = _nThresh;
	nRoots = _nRoots;
	divideByCSFs = _divideByCSFs;
	thresh = new EnergyType[nThresh];
	sum = new EnergyType[nThresh*nRoots];
	countConf = new INT[nThresh];
	countCSF = new INT[nThresh];

	memset(thresh, 0, nThresh*sizeof(EnergyType));	
	memset(sum, 0, nThresh*nRoots*sizeof(EnergyType));
	memset(countConf, 0, nThresh*sizeof(INT));
	memset(countCSF, 0, nThresh*sizeof(INT));
}

PTSum::~PTSum()
{
	delete thresh;
	delete sum;
	delete countConf;
	delete countCSF;
}



ostream& operator<<(ostream & s, const PTSum &ptsum)
{
	s << "      threshold/mH      sel. conf.       sel. CSFs";
	for ( INT i=0 ; i<ptsum.nRoots ; i++ )
		s << "     energy sum/mH";
	s << endl;
	if ( ptsum.nRoots>1 )
	{
		s << "                                                    ";
		for ( INT i=0 ; i<ptsum.nRoots ; i++ )
			s << "        root #" << setw(2) << i+1 << "  ";
		s << endl;
	}
	s << "============================================================";
	for ( INT i=0 ; i<ptsum.nRoots ; i++ )
		s << "==================";
	s << endl;
	for ( INT i=0 ; i<ptsum.nThresh ; i++ )
	{
		s << setw(18) << setiosflags(ios::scientific) << setprecision(6) << ptsum.thresh[i]*1e3
			<< setw(16) << ptsum.countConf[i]
			<< setw(16) << ptsum.countCSF[i];
		for ( INT j=0 ; j<ptsum.nRoots ; j++ )
			s << setw(18) << setiosflags(ios::scientific) << setprecision(6)
				<< ptsum.sum[i + j*ptsum.nThresh]*1e3;
		s << endl;
	}
	s << endl;
	return s;
}
