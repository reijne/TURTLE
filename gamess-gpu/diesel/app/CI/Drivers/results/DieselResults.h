//***********************************************************************
//
//	Name:			DieselResults.h
//
//	Description:	stores input for diesel driver 
//
//	Author:			Michael Hanrath
//
//	Version:		0.0
//
//	Date:			09.11.1998
//
//
//
//
//
//***********************************************************************


#ifndef __DieselResults_H
#define __DieselResults_H


#include "../../../../config.h"


#include "../../../../lib/Container/SLList.h"
#include "../../../../lib/Container/String.h"
#include "../../../../lib/Math/MatrixVector/Matrix.h"
#include "../../../../lib/Math/etc/Number.h"

#include <string>


//FD class ostream;

class DieselResults {
public:
	DieselResults(ostream &s, Number<double> Thresh, INT Root, 
		INT noWave, INT noHeaders, String MRPTThresh);
	~DieselResults();



	friend ostream & operator << (ostream &, const DieselResults &);

protected:
INT	nRoots;
INT	nThresholds;


struct TThreshRootData {
	Number<double>	ECIvar;
	Number<double>	RefCISqr;
	Number<double>	EPTsum_EN;
	Number<double>	EPTsum_ENweighted;
	Number<double>	EPTsum_MRMP2;
	Number<double>	EPTsum_MRMP3;
	Number<double>	overlap;
	INT	rootNr_ref2sel;
	INT	rootNr_sel2ref;
	Number<double>	MRCIextPol[9];
	Number<double>	fullCIextPol[12];
        std::string	remark;
	};

struct TThreshData {
	Number<double>	threshold;
	Number<INT>	nConfs;
	Number<INT> nSAFs;
	TThreshRootData	*root;
	Matrix<double>	overlap;
	Matrix<double>	overlapMax;
	Matrix<double>	overlapRenorm;

	};
	
TThreshData	*threshData;

struct TRootData {
	Number<double>	RefE;
	SLList<String>	waveFunction;
	};

TRootData	*rootData;
String	waveFunctionThreshold;

static const double imp;


Number<double> Thresh;
INT Root;
String Category;
INT	MRMP;
INT noWave;
INT noHeaders;
};




#endif
