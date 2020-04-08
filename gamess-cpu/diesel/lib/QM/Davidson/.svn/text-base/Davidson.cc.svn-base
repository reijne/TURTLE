//***********************************************************************
//
//	Name:			Davidson.cc
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



#include "Davidson.h"
#include "../RepresentationMatrices/CIVectors.h"
#include "../MRCIMatrix/MRCIMatrix.h"

#include "../../Math/MatrixVector/Matrix.h"


#include "../IO/TimeTicks.h"

#include <stdlib.h>
#include <fstream>
#include <iomanip>

#include "../../../config.h"

#include "../../Math/FortranNumeric/FortranEigenProblems.h"
#include "../MRTree/Diag/InternalConfsDiag.h"

#include "../IO/Verbosity.h"

using namespace std;

template <class CIMatrixType, class CIVectorType>
Davidson<CIMatrixType, CIVectorType>::Davidson(INT _dim, INT _RefMatDim, 
	Roots &_roots, const DavidsonMatrix<CIMatrixType, CIVectorType> &_davidsonMatrix, 
	INT restart,
	INT	rootHoming) :
		rootHoming(rootHoming)
{
	dim = _dim;
	RefMatDim = _RefMatDim;
	roots = &_roots; 
INT	mode = DiskBuffer::deleteOnClose;
	if ( !restart )
		mode |= DiskBuffer::truncateOnOpen;

	tBuf = NULL;
	bBuf = new DiskBuffer(dim*sizeof(CIVectorType), "Davidson_b.dat", mode);
	AbBuf = new DiskBuffer(dim*sizeof(CIVectorType), "Davidson_Ab.dat", mode);
	
	Aproj = lambda = alpha = NULL;
	RefMatEigenvectors = NULL;
	
	oldLambda = (CIMatrixType *) malloc(sizeof(CIMatrixType)*roots->getNumberOfRoots());
	memset(oldLambda, 0, sizeof(CIMatrixType)*roots->getNumberOfRoots());

	rootMap = new INT[roots->getNumberOfRoots()];
	for ( INT i=0 ; i<roots->getNumberOfRoots() ; i++ )
		rootMap[i] = roots->getRootNumber(i);
	
	davidsonMatrix = &_davidsonMatrix; 

	diagonal = new CIMatrixType[dim];
	davidsonMatrix->getDiagonal(diagonal);
		
}



template <class CIMatrixType, class CIVectorType>
Davidson<CIMatrixType, CIVectorType>::~Davidson()
{
	delete bBuf;
	delete AbBuf;
	if ( tBuf )
		delete tBuf;
		
	delete[] rootMap;
		
	if ( Aproj ) 
		free(Aproj);
	if ( lambda )
		free(lambda);
	if ( alpha )
		free(alpha);
	if ( RefMatEigenvectors )
		delete[] RefMatEigenvectors;
		

	if ( diagonal )
		delete[] diagonal;
	free(oldLambda);
}



template <class CIMatrixType, class CIVectorType>
void	Davidson<CIMatrixType, CIVectorType>::calcRefMat(
	CIMatrixType *pRefMat, CIMatrixType *lambda, CIMatrixType *alpha)
{
CIMatrixType	*A = new CIMatrixType[RefMatDim*RefMatDim];
	davidsonMatrix->getRefMatrix(A);
	{
	INT	k = 0;
		for ( INT i=0 ; i<RefMatDim ; i++ )
			for ( INT j=0 ; j<=i ; j++ )
				pRefMat[k++] = A[i*RefMatDim + j];
	}
	delete[] A;

	cout << "RefMatDim=" << RefMatDim << endl;
	if ( verbosity.isActive(Verbosity::RefMat) )
	{
		cout << "reference matrix:" << endl;
		for ( INT i=0 ; i<RefMatDim ; i++ )
		{
			for ( INT j=0 ; j<RefMatDim ; j++ )
				if ( j<=i )
					cout << setw(14) << pRefMat[i*(i+1)/2 + j] << " ";
				else
					cout << "\t";
			cout << endl;
		}
		cout << endl;
	}

	RefMatEigenvectors = new CIVectorType[roots->getNumberOfRoots()*RefMatDim];

	hqrii1(
		&RefMatDim, &RefMatDim, pRefMat, lambda, alpha, &RefMatDim
	);

	for ( INT i=0 ; i<roots->getNumberOfRoots() ; i++ )
		for ( INT j=0 ; j<RefMatDim ; j++ )
			RefMatEigenvectors[RefMatDim*i + j] = 
				alpha[RefMatDim*roots->getRootNumber(i)+j];

}


template <class CIMatrixType, class CIVectorType>
void	Davidson<CIMatrixType, CIVectorType>::start()
{
CIMatrixType	*pRefMat = new CIMatrixType[RefMatDim*(RefMatDim+1)/2];
CIMatrixType	*lambda = new CIMatrixType[RefMatDim];
CIMatrixType	*alpha = new CIMatrixType[RefMatDim*RefMatDim];

	calcRefMat(pRefMat, lambda, alpha);


	if ( verbosity.isActive(Verbosity::RefMatEigenValues) )
	{
		cout << "eigenvalues of reference matrix:" << endl;
		for ( INT i=0 ; i<RefMatDim ; i++ )
			cout << lambda[i] << endl;	
		cout << endl;
	}


	if ( verbosity.isActive(Verbosity::RefMatEigenVectors) )
	{
		cout << "eigenvectors of reference matrix (columnwise):" << endl;
		for ( INT i=0 ; i<RefMatDim ; i++ )
		{
			for ( INT j=0 ; j<RefMatDim ; j++ )
				cout << setw(14) << alpha[j*RefMatDim + i] << " ";
			cout << endl;
		}
		cout << endl;
	}



	if ( roots->getNumberOfRoots()>RefMatDim )
	{
		cout << "too many (" << roots->getNumberOfRoots() << ") roots ordered" << endl;
		exit(1);
	}

//	add some noise to start vectors to insure start of davidson algorithm
//	in pathologic cases
	cout << "checking condition of starting vectors:" << endl;
	for ( INT i=0 ; i<RefMatDim ; i++ )
	{
	double	vecSum = 0;

		for ( INT j=0 ; j<RefMatDim ; j++ )
			vecSum += alpha[i*RefMatDim + j];
		cout << "checking vector #" << i << "...";

		srand(1000);
		if ( fabs(fabs(vecSum)-1)<0.01 )
		{
			lambda[i] -= 0.01*(1.0*rand())/RAND_MAX;

			for ( INT j=0 ; j<RefMatDim ; j++ )
				alpha[i*RefMatDim + j] -= 0.01*(1.0*rand())/RAND_MAX;
			cout << "adding noise." << endl;
		}
		else
			cout << "OK." << endl;
	}
	cout << endl;


INT	mode = DiskBuffer::deleteOnClose | DiskBuffer::truncateOnOpen;
	tBuf = new DiskBuffer(dim*sizeof(CIVectorType), "Davidson_t.dat", mode);
CIVectorType	*BasisVector = new CIVectorType[dim];

	for ( INT j=0 ; j<roots->getNumberOfRoots() ; j++ )
	{
		memset(BasisVector, 0, dim*sizeof(CIVectorType));
		for ( INT i=0 ; i<RefMatDim ; i++ )
			BasisVector[davidsonMatrix->getRefIndex(i)]
				= alpha[roots->getRootNumber(j)*RefMatDim + i];

		roots->setRefRoot(j, lambda[roots->getRootNumber(j)]);

	lambda[roots->getRootNumber(j)] = roots->getRefRoot(j);

		tBuf->put(j, BasisVector);
	}
	delete[] pRefMat;
	delete[] lambda;
	delete[] alpha;

	delete[] BasisVector;

	start(tBuf);

	if ( !rootHoming )
	{
		delete tBuf;
		tBuf = NULL;
	}
}

template <class CIMatrixType, class CIVectorType>
const CIVectorType *	Davidson<CIMatrixType, CIVectorType>::getRefMatEigenvectorP(INT root) const
{
	return RefMatEigenvectors + RefMatDim*root;
}







template <class CIMatrixType, class CIVectorType>
void	Davidson<CIMatrixType, CIVectorType>::start(DiskBuffer *startVectors)
{
	
CIMatrixType	*pRefMat = new CIMatrixType[RefMatDim*(RefMatDim+1)/2];
CIMatrixType	*lambda = new CIMatrixType[RefMatDim];
CIMatrixType	*alpha = new CIMatrixType[RefMatDim*RefMatDim];

	calcRefMat(pRefMat, lambda, alpha);

	delete[] pRefMat;
	delete[] lambda;
	delete[] alpha;


	basisDim = startVectors->getNumberOfObjects();

//	printf("basisDim= %d\n", basisDim);
	cout << "basisDim=" << basisDim << endl;

	bBuf->clear();

Vector<CIVectorType>	*v1;

	v1 = new Vector<CIVectorType>(dim);

	if ( rootHoming && !tBuf )
	{
	INT	mode = DiskBuffer::deleteOnClose | DiskBuffer::truncateOnOpen;
		tBuf = new DiskBuffer(dim*sizeof(CIVectorType), "Davidson_t.dat", mode);
		for ( INT i=0 ; i<startVectors->getNumberOfObjects() ; i++ )
		{
			startVectors->get(i, v1->getP());
			tBuf->put(i, v1->getP());
		}
	}

	for ( INT i=0 ; i<basisDim ; i++ )
	{	startVectors->get(i, v1->getP());
//		for ( INT j=0 ; j<10 ; j++ )
//			printf("%lf\n", (*v1)[j]);
		orthogonalize(bBuf, v1, i);
	double	Norm = v1->getNorm2();
		*v1 /= Norm;
//		cout << "Norm=" << Norm << endl;
//		for ( INT j=0 ; j<10 ; j++ )
//			printf("%lf\n", (*v1)[j]);
		bBuf->put(i, v1->getP());
	}
//	cout << *v1 << endl;
	delete v1;


	davidsonMatrix->mult(bBuf, AbBuf, 0, basisDim-1);


	calcProjected();

//	printOrtho(bBuf);
}


template <class CIMatrixType, class CIVectorType>
void	Davidson<CIMatrixType, CIVectorType>::restart()
{
CIMatrixType	*pRefMat = new CIMatrixType[RefMatDim*(RefMatDim+1)/2];
CIMatrixType	*lambda = new CIMatrixType[RefMatDim];
CIMatrixType	*alpha = new CIMatrixType[RefMatDim*RefMatDim];

	calcRefMat(pRefMat, lambda, alpha);

	delete[] pRefMat;
	delete[] lambda;
	delete[] alpha;


	if ( bBuf->getNumberOfObjects()>AbBuf->getNumberOfObjects() )
		bBuf->deleteLast();
	basisDim = bBuf->getNumberOfObjects();
//	printf("basisDim= %d\n", basisDim);
	cout<<"basisDim= "<<basisDim<<endl;

	calcProjected();
}


template <class CIMatrixType, class CIVectorType>
void	Davidson<CIMatrixType, CIVectorType>::calcProjected(INT start)
{
TimeTicks	t;
	t.start();

	cout << "calculating projected matrix..." << flush;

	allocateProjected(basisDim);

INT	nBuf = roots->getNumberOfRoots();
Vector<CIVectorType>       *v1 = new Vector<CIVectorType>[nBuf];
Vector<CIVectorType>       *v2 = new Vector<CIVectorType>[nBuf];
for ( INT i=0; i<nBuf; i++)
{
	v1[i].setDim(dim);
	v2[i].setDim(dim);
}

INT	blocks = 1 + (basisDim-1) / nBuf;

//	cout << "basisDim=" << basisDim << endl;


	for ( INT iblock=start/nBuf ; iblock<blocks ; iblock++ )
	{
		for ( INT i=0 ; i<nBuf && iblock*nBuf+i<basisDim ; i++ )
			bBuf->get(iblock*nBuf + i, v1[i].getP());
	
		for ( INT jblock=0 ; jblock<=iblock ; jblock++ )
		{
			for ( INT i=0 ; i<nBuf && jblock*nBuf+i<basisDim ; i++ )
				AbBuf->get(jblock*nBuf + i, v2[i].getP());

			for ( INT i=0 ; i<nBuf && iblock*nBuf+i<basisDim ; i++ )
			{
			INT	ii = iblock*nBuf + i;
				for ( INT j=0 ; j<=(iblock==jblock ? i : nBuf-1) && jblock*nBuf+j<basisDim ; j++ )
				{
				INT	jj = jblock*nBuf + j;
				
				
//					cout << "ii=" << ii << ", jj=" << jj << ": " << v1[i]*v2[j] << endl;
					Aproj[ii*(ii+1)/2+jj] = v1[i]*v2[j];
				
				}
			}



/*
		INT	is = start;

			for ( INT i=iStart ; i<basisDim ; i++ )
			{
				for ( INT j=0 ; j<=i ; j++ )
				{
					if ( i>=iblock*nBuf && i<(iblock+1)*nBuf &&
						 j>=jblock*nBuf && j<(jblock+1)*nBuf )
					{
						cout << "i=" << i << ", j=" << j << endl;
						AbBuf->get(j, v2->getP());
		//				Aproj[j*basisDim + i] = Aproj[i*basisDim + j] = *v1 * *v2;
			TimeTicks	t;
						t.start();
						Aproj[is++] = *v1 * *v2;
						t.stop();
						cout << t << endl;
		//				cout << *v1 * *v2 << endl;
					}
					else
						is++;
				}
			}
			
*/
		}
	}


	delete[] v1;
	delete[] v2;

	t.stop();
	cout << t << endl;
}


template <class CIMatrixType, class CIVectorType>
void	Davidson<CIMatrixType, CIVectorType>::iterate(
	INT maxIters,
	ConvergenceCriterion crit,
	double precision)
{
INT	nconv;
	for ( INT k=0 ; k<roots->getNumberOfRoots() ; k++ )
		roots->setConvergenceStatus(k, Roots::NotConverged);

	cout << setw(5) << 0;	
	for ( INT j=0 ; j<roots->getNumberOfRoots() ; j++ )
		cout << setw(22) << "---" <<
			(roots->getConvergenceStatus(j)==Roots::Converged) ? " (*) " : "     ";
	cout << endl;	

	for ( INT i=0 ; i<roots->getNumberOfRoots() ; i++ )
		lambda[i] = roots->getRefRoot(i);

INT	status = 0;
	for ( INT i=0 ; i<maxIters ; i++ )
	{
		if ( doIteration(crit, precision, nconv, i) )
		{
			status = -1;
			break;
		}


		cout << setw(5) << i+1;	
		for ( INT j=0 ; j<roots->getNumberOfRoots() ; j++ )
			cout << setw(22) << setprecision(11) << lambda[rootMap[j]] <<
				((roots->getConvergenceStatus(j)==Roots::Converged) ? " (*) " : "     ");
		cout << endl;	


//		cout << "iteration " << i+1 << " completed" << endl;
//		cout << nconv << " out of ";
//		cout << roots->getNumberOfRoots() << " roots converged" << endl;
		if ( nconv==roots->getNumberOfRoots() )
		{
			status = 1;
			break;
		}
	}

	for ( INT k=0 ; k<roots->getNumberOfRoots() ; k++ )
		roots->setRoot(k, lambda[rootMap[k]]);

	cout << "*************************************************" << endl;
	cout << "***                                           ***" << endl;

	switch ( status ) {
	case -1:
		cout << "***       convergence was not achieved:       ***" << endl;
		cout << "***         correction vector colinear        ***" << endl;
		break;
	
	case 0:
		cout << "***       convergence was not achieved:       ***" << endl;
		cout << "***           too many iterations             ***" << endl;
		break;

	case 1:
		cout << "***           all roots converged             ***" << endl;
		break;
	}
	cout << "***                                           ***" << endl;
	cout << "*************************************************" << endl << endl;




	if ( diagonal )
		delete[] diagonal;
	diagonal = NULL;

//	if ( status==1 )
	{
	//	calculate eigenvectors
Vector<CIVectorType>       *ev = new Vector<CIVectorType>[roots->getNumberOfRoots()];
for ( INT i=0; i<roots->getNumberOfRoots(); i++)
{
	ev[i].setDim(dim);
}
Vector<CIVectorType>	*bv = new Vector<CIVectorType>(dim);
DiskBuffer	*evBuf = new DiskBuffer(dim*sizeof(CIVectorType), "Eigenvectors.dat",
							DiskBuffer::truncateOnOpen | DiskBuffer::noTempDir);

		for ( INT i=0 ; i<roots->getNumberOfRoots() ; i++ )
			ev[i].clear();

		for ( INT j=0 ; j<basisDim ; j++ )
		{
			bBuf->get(j, bv->getP());
			for ( INT i=0 ; i<roots->getNumberOfRoots() ; i++ )
				ev[i] += *bv * alpha[rootMap[i]*basisDim + j];
		}
		for ( INT i=0 ; i<roots->getNumberOfRoots() ; i++ )
			evBuf->put(i, ev[i].getP());

		delete evBuf;
		delete[] ev;
		delete bv;
	}

}




template <class CIMatrixType, class CIVectorType>
INT	Davidson<CIMatrixType, CIVectorType>::doIteration(
	ConvergenceCriterion crit,
	double precision,
	INT &nconv, INT iter)
{

/*	cout << "projected matrix:" << endl;
	for ( INT i=0 ; i<basisDim ; i++ )
	{
		for ( INT j=0 ; j<basisDim ; j++ )
			if ( j<=i )
				cout << Aproj[i*(i+1)/2 + j] << " ";
			else
				cout << "\t";
		cout << endl;
	}
	cout << endl;
*/

Vector<CIVectorType>	*x;
Vector<CIVectorType>	*y = new Vector<CIVectorType>(dim);


INT	bInc = 0;
DiskBuffer	*yBuf = new DiskBuffer(dim*sizeof(CIVectorType), 
	"Davidson_y.dat");


//*****************************************************************************
//
//                               -->
//  calculate correction vectors  f
//                                 k
//
//                                    L
//                                  ------
//   (-->)             1             \       ( -->  )                   (-->)
//   ( f )   = -----------------   *  >      (alpha )   (A - lambda ) * ( b )
//   (  k)      lambda  - (A)        /       (     k)              k    (  i)
//        j          k      jj      ------           i                       j 
//                                    i=1
//
//  with
//
//  A      : Hamilton Matrix
//
//  lambda : k-th eigenvalue of projected matrix
//        k
//
//   -->
//  alpha  : k-th eigenvector of projected matrix
//       k
//
//  -->
//   b     : basis vector of projected matrix
//    i
//
//  M      : number of desired roots          = "roots"
//  L      : dimension of projected matrix    = "basisDim"
//  N      : dimension of A (full ci-matrix)  = "dim"
//
//  k = {1, ..., M}
//  i = {1, ..., L}
//  j = {1, ..., N}
//
//
//
//   (-->) 
//	 ( x )  denotes the k-th component of vector x
//        k
//
//*****************************************************************************


	nconv = 0;

	/////////////////////////////////////////////////////////////////////
	//	optimized for minimum memory usage
	//

	{	
	CIMatrixType	*A = new CIMatrixType[basisDim*(basisDim+1)/2];

		memcpy(A, Aproj, basisDim*(basisDim+1)/2*sizeof(CIMatrixType));


		hqrii1(
			&basisDim, &basisDim, A, lambda, alpha, &basisDim
		);

		delete[] A;


		if ( rootHoming )
		{
		const INT	nTest = tBuf->getNumberOfObjects();
		Matrix<CIVectorType>	U(basisDim, nTest);
			{
			Vector<CIVectorType>	b(dim);
			Vector<CIVectorType>	t(dim);


				for ( INT k=0 ; k<basisDim ; k++ )
				{
					bBuf->get(k, b.getP());
					for ( INT j=0 ; j<nTest ; j++ )
					{
						tBuf->get(j, t.getP());
						U(k, j) = b*t;
					}
				}
			}
			Matrix<CIMatrixType>	alphaM(basisDim, basisDim, alpha);
			Matrix<CIMatrixType>	M(basisDim, nTest);
				M = alphaM*U;
//			cout << "alpha=" << alphaM << endl;
//			cout << "U=" << U << endl;
			cout << "M=" << M << endl;
			
			
			// search for maximum overlap
 			INT* used = new INT[basisDim];
			memset(used, 0, basisDim*sizeof(INT));
			for ( INT j=0 ; j<nTest ; j++ )
			{
			CIMatrixType	max = 0;
			INT	argmax = 0;
				for ( INT k=0 ; k<basisDim ; k++ )
				{
					if ( fabs(M(k, j))>max && !used[k] )
					{
						max = fabs(M(k, j));
						argmax = k;
					}
				}
				used[argmax] = 1;
				cout << j << ": argmax=" << argmax << endl;
				rootMap[j] = argmax;
			}
			delete used;
		}


//	add some noise to start vectors to insure start of davidson algorithm
//	in pathologic cases
		if ( !iter )
		{
			if ( RefMatDim==1 )
			{
				for ( INT i=0 ; i<RefMatDim ; i++ )
						lambda[i] -= 0.01*(1.0*rand())/RAND_MAX;

				for ( INT i=0 ; i<RefMatDim ; i++ )
				{
					for ( INT j=0 ; j<RefMatDim ; j++ )
						alpha[j*RefMatDim + i] -= 0.01*(1.0*rand())/RAND_MAX;
				}
			}
		}
	}

INT	kk = 0;
	for ( INT k=0 ; k<roots->getNumberOfRoots() ; k++ )
	{


		if ( roots->getConvergenceStatus(k)==Roots::Converged )
		{
			nconv++;
			continue;
		}

		if ( crit == Energy )
			if ( fabs(lambda[rootMap[k]] -
					oldLambda[k])<precision )
			{
				roots->setConvergenceStatus(k, Roots::Converged);
				nconv++;
				continue;
			}


		x = new Vector<CIVectorType>(dim);
	
	CIMatrixType	lambdaK = lambda[rootMap[k]];
	
		
		
		y->clear();
		
		for ( INT i=0 ; i<basisDim ; i++ )
		{
		CIMatrixType	alphaKI = alpha[rootMap[k]*basisDim + i];

			{
			// calculate alpha_k*A*b_i
				AbBuf->get(i, x->getP());
		
			CIMatrixType	*pADiag = diagonal;
			CIVectorType	*pAb = x->getP();
			CIVectorType	*pf = y->getP();

//		cout << "lambdaK=" << lambdaK << endl;
//		cout << "*pADiag=" << *pADiag << endl;
//		cout << "alphaKI=" << alphaKI << endl;

				for ( INT j=0 ; j<dim ; j++ )
					*pf++ += alphaKI/(lambdaK - *pADiag++) * *pAb++;
			}
			{
			// calculate alpha_k*lambda_k*b_i
				bBuf->get(i, x->getP());

			CIMatrixType	*pADiag = diagonal;
			CIVectorType	*pb = x->getP();
			CIVectorType	*pf = y->getP();

				for ( INT j=0 ; j<dim ; j++ )
					*pf++ -= alphaKI/(lambdaK - *pADiag++) * lambdaK * *pb++;
			}
		}
		delete x;
		if ( crit == CorrectionVector )
			if ( fabs(y->getNorm2())<precision )
			{
				roots->setConvergenceStatus(k, Roots::Converged);
				nconv++;
				continue;
			}

		oldLambda[k] = lambda[rootMap[k]];

	
		y->normalize2();
		yBuf->put(kk++, y->getP());
	}
	
	
//	orthogonalize new vectors
	orthogonalize(bBuf, yBuf);
		
	for ( INT k=0 ; k<yBuf->getNumberOfObjects() ; k++ )
	{
		yBuf->get(k, y->getP());

//		printOrtho(bBuf);

//	check if new vector is worth to be appended to the basis vectors
//	cout << *y << endl;
	double	t = y->getNorm2();
//	cout << "Norm=" << t << endl;
		if ( t>SelThreshold )
			bBuf->put(basisDim + bInc++, y->getP());
	}
//	printOrtho(bBuf);
	delete y;
	delete yBuf;


	if ( !bInc )
		return nconv<roots->getNumberOfRoots();
	
	
//	calculate A*b
	davidsonMatrix->mult(bBuf, AbBuf, basisDim, basisDim+bInc-1);

	basisDim += bInc;

//	cout << bInc << " vectors appended to basis." << endl;
//	cout << "The basis now consists of " << basisDim << " vectors." << endl;



	allocateProjected(basisDim);

	calcProjected(basisDim-bInc);

	return 0;
}	





template <class CIMatrixType, class CIVectorType>
void	Davidson<CIMatrixType, CIVectorType>::allocateProjected(INT dim)
{
	Aproj = (CIMatrixType *) realloc(Aproj, sizeof(CIMatrixType)*dim*(dim+1)/2);
	lambda = (CIMatrixType *) realloc(lambda, sizeof(CIMatrixType)*dim);
	alpha = (CIMatrixType *) realloc(alpha, sizeof(CIMatrixType)*dim*dim);
}


template <class CIMatrixType, class CIVectorType>
void	Davidson<CIMatrixType, CIVectorType>::orthogonalize(DiskBuffer *basisBuf, DiskBuffer *yBuf)
{
TimeTicks	t;
	t.start();

	cout << "orthogonalizing new vectors..." << flush;


Vector<CIVectorType>	base(dim);

Vector<CIVectorType>        *y = new Vector<CIVectorType>[yBuf->getNumberOfObjects()];
for ( INT i=0; i<yBuf->getNumberOfObjects(); i++)
{
	y[i].setDim(dim);
}

	for ( INT i=0 ; i<yBuf->getNumberOfObjects() ; i++ )
		yBuf->get(i, y[i].getP());

	// orthogonalize on basis
	for ( INT i=0 ; i<basisBuf->getNumberOfObjects() ; i++ )
	{
		basisBuf->get(i, base.getP());

//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
//	caution:
//	numerical inaccuracies may cause problems
		for ( INT j=0 ; j<yBuf->getNumberOfObjects() ; j++ )
		{
		CIVectorType	SP = base * y[j];
		CIVectorType	*p1 = y[j].getP();
		CIVectorType	*p2 = base.getP();
			for ( INT k=0 ; k<dim ; k++ )
				*p1++ -= SP * *p2++;
		}
	}

INT*	toSmall = new INT[yBuf->getNumberOfObjects()];
	memset(toSmall, 0, yBuf->getNumberOfObjects()*sizeof(INT));
	
	// orthogonalize on each other
	for ( INT i=0 ; i<yBuf->getNumberOfObjects() ; i++ )
	{
	double	norm2 = y[i].getNorm2();

		if ( norm2>SelThreshold )
			y[i] /= norm2;
		else
			toSmall[i] = 1;
			
		if ( toSmall[i] )
			continue;
		for ( INT j=i+1 ; j<yBuf->getNumberOfObjects() ; j++ )
		{
			if ( toSmall[j] )
				continue;

		CIVectorType	SP = y[i] * y[j];
		CIVectorType	*p1 = y[j].getP();
		CIVectorType	*p2 = y[i].getP();
			for ( INT k=0 ; k<dim ; k++ )
				*p1++ -= SP * *p2++;

		}
	}


	for ( INT i=0 ; i<yBuf->getNumberOfObjects() ; i++ )
		yBuf->put(i, y[i].getP());
	
	delete[] y;
	t.stop();
	cout << t << endl;
	delete[] toSmall;
}

template <class CIMatrixType, class CIVectorType>
void	Davidson<CIMatrixType, CIVectorType>::orthogonalize(DiskBuffer *Buf, Vector<CIVectorType> *b, INT n)
{
TimeTicks	t;
	t.start();

	cout << "orthogonalizing new vector..." << flush;


Vector<CIVectorType>	v(dim);

	for ( INT i=0 ; i<n ; i++ )
	{
		Buf->get(i, v.getP());
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
//	caution:
//	numerical inaccuracies may cause problems
	CIVectorType	SP = v * *b;
	CIVectorType	*p1 = b->getP();
	CIVectorType	*p2 = v.getP();
		for ( INT j=0 ; j<dim ; j++ )
			*p1++ -= SP * *p2++;
	}
	t.stop();
	cout << t << endl;
}


template <class CIMatrixType, class CIVectorType>
void	Davidson<CIMatrixType, CIVectorType>::orthogonalize(DiskBuffer *Buf, Vector<CIVectorType> *b)
{
	orthogonalize(Buf, b, Buf->getNumberOfObjects());
}


template <class CIMatrixType, class CIVectorType>
void	Davidson<CIMatrixType, CIVectorType>::printOrtho(DiskBuffer *Buf)
{
INT	n = Buf->getNumberOfObjects();
INT dim = Buf->getObjectSize() / sizeof(CIVectorType);
Vector<CIVectorType>	v1(dim);
Vector<CIVectorType>	v2(dim);

	cout << "orthogonalization matrix:" << endl;
	for ( INT i=0 ; i<n ; i++ )
	{
		Buf->get(i, v1.getP());
		for ( INT j=0 ; j<=i ; j++ )
		{
			Buf->get(j, v2.getP());
			cout << v1*v2 << " ";
		}
		cout << endl;
	}
	cout << endl;
}


template class Davidson<float, float>;
//template class Davidson<float, double>;
//template class Davidson<double, float>;
template class Davidson<double, double>;

template <class CIMatrixType, class CIVectorType>
const double Davidson<CIMatrixType, CIVectorType>::SelThreshold = 1e-3;
