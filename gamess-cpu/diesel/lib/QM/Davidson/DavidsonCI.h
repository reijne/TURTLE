//***********************************************************************
//
//	Name:			DavidsonCI.h
//
//	Description:	
//
//	Author:			Michael Hanrath
//
//	Version:		0.0
//
//	Date:			14.08.1998
//
//
//
//
//***********************************************************************


#ifndef __DavidsonCI_H
#define __DavidsonCI_H


#include "../../../config.h"

#include "Davidson.h"

#include "../Configuration/Configuration.h"

class NExternalsDiag;

template <class CIMatrixType, class CIVectorType> class DavidsonMatrix;

template <class CIMatrixType, class CIVectorType>
class DavidsonCI : public Davidson<CIMatrixType, CIVectorType> {
public:
enum IterationMode	{ CI, ACPF, AQCC };
	DavidsonCI(
		INT dim, 
		INT RefMatDim, 
		Roots &roots, 
		const NExternalsDiag *,
		const DavidsonMatrix<CIMatrixType, CIVectorType> &,
		INT restart,
		INT	rootHoming, 
		IterationMode method = CI);


	~DavidsonCI();

	//	initialize basis with eigenvectors from reference space
	void	start();

	//	initialize basis with eigenvectors from previous calculation
	void	start(
		const char *ConfTreeFileName,
		const char *EigenvectorsFileName);

	//	initialize basis with given vectors
	void	start(DiskBuffer *startVectors);

	//	restart from "Davidson_b.dat" and "Davidson_Ab.dat"
	void	restart();
        void    iterate(
                INT maxIters = 20,
                ConvergenceCriterion crit = Energy,
                double precision = 1e-5);


private:
	INT	doIteration(
		typename DavidsonCI<CIMatrixType,CIVectorType>::ConvergenceCriterion crit,
		double precision,
		INT &nconv, INT iter);

	INT	doIterationSzalay(
		typename DavidsonCI<CIMatrixType,CIVectorType>::ConvergenceCriterion crit,
		double precision,
		INT &nconv, INT iter);


	void	calcExtPart(INT dim, CIVectorType *p);
	INT		checkInactiveHole(Configuration<MOType> &, const MRMOs *);



const NExternalsDiag 		*extDiag;	//	needed for Hamilton-Matrix multiplication

IterationMode	method;			//	iteration method: (CI, ACPF, AQCC) 


CIVectorType	*ci2Ext;		//	sum over CSF coefficients of external
								//	configurations and such with an hole
								//	in an MO occupied in every reference
								//	(an inactive one)
								
INT				internalCSFs;	//	number of CSFs in internal space
INT				*inactiveHole;	//	flag if configuration has a hole in
								//	an inactive MO
};





#endif
