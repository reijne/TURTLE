//***********************************************************************
//
//	Name:			PTSum.h
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

#ifndef __PTSum_h
#define __PTSum_h

#include "../../../../config.h"

#include "../EnergyType.h"

#include <iostream>
using std::ostream;

class PTSum {
public:
	PTSum(EnergyType SelectionThreshold, INT divideByCSFs, INT nRoots = 1);
	PTSum(INT nThresh, const EnergyType *e, INT divideByCSFs, INT nRoots = 1);
	PTSum(INT nThresh, INT divideByCSFs, INT nRoots = 1);
	~PTSum();


	void	count(EnergyType e, INT CSFs, INT selectAll);
	void	count(const EnergyType *e, INT CSFs, INT selectAll);
	
	void	setE(INT threshNr, INT rootNr, EnergyType e);
	void	setConfSAF(INT threshNr, EnergyType _thresh, INT Confs, INT SAFs);

//--------------------------------------------------------------------------

	friend ostream& operator<<(ostream & s, const PTSum &);

//--------------------------------------------------------------------------

private:
INT	nThresh;			// number of thresholds
INT	nRoots;				// number of roots
EnergyType	*thresh;	// thresholds
EnergyType	*sum;		// energy sum resulting from not selected configurations
INT	*countConf;			// number of configurations selected on threshold
INT	*countCSF;			// number CSFs selected on threshold
INT	divideByCSFs;		// flag if energy estimation are to be divided by #CSFs

};


#include <math.h>


inline
void	PTSum::count(EnergyType e, INT CSFs, INT selectAll)
{
	for ( INT i=0 ; i<nThresh ; i++ )
		if ( selectAll || fabs(e)/CSFs>=thresh[i] )
		{
			countConf[i]++;
			countCSF[i] += CSFs;
		}
		else
			sum[i] += e;
}


inline
void	PTSum::count(const EnergyType *e, INT CSFs, INT selectAll)
{

	for ( INT i=0 ; i<nThresh ; i++ )
	{
	INT	sel = selectAll;
		if ( !sel )
			for ( INT j=0 ; j<nRoots ; j++ )
//				if ( fabs(e[j])>=thresh[i] )
				if ( fabs(e[j])/CSFs>=thresh[i] )
				{
					sel = 1;
					break;
				}
		if ( sel )
		{
			countConf[i]++;
			countCSF[i] += CSFs;
		}
		else
			for ( INT j=0 ; j<nRoots ; j++ )
				sum[i + j*nThresh] += e[j];
	}
}


inline
void	PTSum::setE(INT threshNr, INT rootNr, EnergyType e)
{
	sum[threshNr + rootNr*nThresh] = e;
}

inline
void	PTSum::setConfSAF(INT threshNr, EnergyType _thresh, INT Confs, INT SAFs)
{
	thresh[threshNr] = _thresh;
	countConf[threshNr] = Confs;
	countCSF[threshNr] = SAFs;
}



#endif
