//***********************************************************************
//
//	Name:			Davidson.h
//
//	Description:	implements the Davidson-Liu Algorithm for
//					calculation of the first few eigenvalues and
//					eigenvectors of VERY large (dim>1e6) matrices
//
//	Author:			Michael Hanrath
//
//	Version:		0.0
//
//	Date:			04.12.1996
//
//	Ref:		1)	E. R. Davidson:
//					Journal of Computational Physics (1975), 17, 87
//
//				2)	B. Liu
//					Numerical Algorithms in Chemistry:
//					Algebraic Methods,
//					LBL-8158 Lawrence Berkeley Laboratory
//					Editors: C. Moler and I. Shavitt
//					
//
//
//
//***********************************************************************


#ifndef __DAVIDSON_H
#define __DAVIDSON_H

#include "../../../config.h"

#include "../../Container/DiskBuffer.h"
#include "Roots.h"
#include "../../Math/MatrixVector/Vector.h"
#include "../RepresentationMatrices/CIVectors.h"


#include "DavidsonMatrix.h"

template <class CIMatrixType, class CIVectorType> class DavidsonMatrix;

template <class CIMatrixType, class CIVectorType>
class Davidson {
public:
enum Mode { direct, storeMatrix };

	Davidson(
		INT dim, 
		INT RefMatDim, 
		Roots &roots, 
		const DavidsonMatrix<CIMatrixType, CIVectorType> &, 
		INT restart,
		INT	rootHoming);
				
	virtual ~Davidson();
	
enum ConvergenceCriterion { Energy, ProjectedEigenvector, CorrectionVector };
	
	//	initialize basis with eigenvectors from reference space
	void	start();

	
	//	initialize basis with given vectors
	void	start(DiskBuffer *startVectors);
	
	//	restart from "Davidson_b.dat" and "Davidson_Ab.dat"
	void	restart();
	
	//	iterate until convergence reached or maximum number of
	//	iterations exceeded
	void	iterate(
		INT maxIters = 20,
		ConvergenceCriterion crit = Energy,
		double precision = 1e-5);



	const CIVectorType *	getRefMatEigenvectorP(INT root) const;

	
static const double	SelThreshold;

protected:
	void	allocateProjected(INT dim);
	void	calcProjected(INT start = 0);
	void	orthogonalize(DiskBuffer *basisBuf, DiskBuffer *yBuf);
	void	orthogonalize(DiskBuffer *bBuf, Vector<CIVectorType> *b);
	void	orthogonalize(DiskBuffer *bBuf, Vector<CIVectorType> *b, INT n);
	void	printOrtho(DiskBuffer *Buf);
	INT		doIteration(
			ConvergenceCriterion crit,
			double precision,
			INT &nconv, INT iter);
	void	calcRefMat(CIMatrixType *pRefMat, CIMatrixType *lambda, CIMatrixType *alpha);
	
	
const DavidsonMatrix<CIMatrixType, CIVectorType>
				*davidsonMatrix;	//	matrix to be diagonalized

INT				rootHoming;			//	flag if root homing is active


INT				dim;			//	dimension of ci-matrix

INT 			RefMatDim;		//	dimension of reference space ci-matrix

Roots			*roots;			//	information about desired roots

INT				*rootMap;		//	root mapping from ordered number --> real number

CIMatrixType	*diagonal;		//	diagonal of A

INT				basisDim;		//	number of basis vectors

DiskBuffer		*tBuf;			//	test vector buffer (for root homing)

DiskBuffer		*bBuf;			//	basis vector buffer

DiskBuffer		*AbBuf;			//	(A * basis vector) buffer

CIMatrixType	*Aproj;			//	projected matrix A
								//	(symmetric, stored in lower triangular)

CIMatrixType	*lambda;		//	eigenvalues of projected Matrix A
CIMatrixType	*oldLambda;		//	previous eigenvalues of projected Matrix A

CIMatrixType	*alpha;			//	eigenvectors of projected Matrix A


CIVectorType	*RefMatEigenvectors; // reference matrix eigenvectors

};


#endif
