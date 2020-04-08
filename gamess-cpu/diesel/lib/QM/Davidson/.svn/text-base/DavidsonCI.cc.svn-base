//***********************************************************************
//
//	Name:			DavidsonCI.cc
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


#include "DavidsonCI.h"
#include "../RepresentationMatrices/CIVectors.h"
#include "../MRCIMatrix/MRCIMatrix.h"

#include "../IO/TimeTicks.h"

#include <stdlib.h>
#include <fstream>
#include <iomanip>

#include "../MO/MOMapping.h"

#include "../MRTree/Sel/NExternalsSel.h"
#include "../MRTree/Diag/NExternalsDiag.h"
#include "../MRTree/Diag/InternalConfsDiag.h"
#include "../MRTree/MRTreeIterator.h"
#include "../Configuration/ConfigurationSAFNr.h"
#include "../MRTree/Sel/extMOsSel.h"

#include "../../Math/FortranNumeric/FortranEigenProblems.h"

#include "../IO/Verbosity.h"

using namespace std;

template <class CIMatrixType, class CIVectorType>
DavidsonCI<CIMatrixType, CIVectorType>::DavidsonCI(
	INT _dim,
	INT _RefMatDim, 
	Roots &_roots,
	const NExternalsDiag *_extDiag,
	const DavidsonMatrix<CIMatrixType, CIVectorType> &_davisonMatrix,
	INT restart,
	INT	rootHoming,
	IterationMode _method) :
	Davidson<CIMatrixType, CIVectorType>( _dim, _RefMatDim, _roots, _davisonMatrix, restart, rootHoming)
{
	ci2Ext = NULL;
	method = _method;


	extDiag = _extDiag;

/////////////////////////////////////////////////////////////////////////////
//*
//*
	//	number of internal CSFs
	internalCSFs = 0;
	{
	const InternalConfsDiag	*diag = (*extDiag)[0];
	
		for ( INT i=0 ; i<diag->getNumberOfElements() ; i++ )
		{
		Configuration<MOType>	conf = *((*diag)[i]);
			internalCSFs += extDiag->getNumberOfSpinAdaptedFunctions(
				conf.getNumberOfOpenShells());
		}
	}


	inactiveHole = new INT[internalCSFs];
	memset(inactiveHole, 0, internalCSFs*sizeof(INT));
	{
	INT	icsf = 0;
	INT	csfs = 0;
	const InternalConfsDiag	*diag = (*extDiag)[0];
	
		for ( INT i=0 ; i<diag->getNumberOfElements() ; i++ )
		{
		Configuration<MOType>	conf = *((*diag)[i]);
			csfs = extDiag->getNumberOfSpinAdaptedFunctions(
				conf.getNumberOfOpenShells());
				
			
//			cout << conf << ": ";
			if ( checkInactiveHole(conf, extDiag->getMRMOs()) )
			{
				for ( INT j=0 ; j<csfs ; j++ )
					inactiveHole[icsf++] = 1;
//				cout << "yes";
			}
			else
				icsf += csfs;
//			cout << endl;
		}
	}

//*
//*
/////////////////////////////////////////////////////////////////////////////

/*
	diagonal = new CIMatrixType[dim];
	cout << "vor diag" << endl;
	mrciMatrix->getDiagonal(diagonal);
	cout << "nach diag" << endl;
*/	
	
/*	cout << "Diagonale: " << endl;
	for ( INT i=0 ; i<RefMatDim ; i++ )
		cout << diagonal[i] << endl;
*/
}


template <class CIMatrixType, class CIVectorType>
DavidsonCI<CIMatrixType, CIVectorType>::~DavidsonCI()
{
	delete[] inactiveHole;

	if ( ci2Ext )
		free(ci2Ext);
}


template <class CIMatrixType, class CIVectorType>
void	DavidsonCI<CIMatrixType, CIVectorType>::start()
{
	Davidson<CIMatrixType, CIVectorType>::start();
}


template <class CIMatrixType, class CIVectorType>
void	DavidsonCI<CIMatrixType, CIVectorType>::start(
	const char *ConfTreeFileName,
	const char *EigenvectorsFileName)
{
DiskBuffer startVectors(this->dim*sizeof(CIVectorType), "Davidson_Basis.dat");

CIVectorType	*BasisVector = new CIVectorType[this->dim];

TimeTicks	ticks;

	cout << "reading configuration tree from previous calculation..." << flush;
	
ifstream	ConfIn(ConfTreeFileName);

MOMapping	moMapping(ConfIn);
NExternalsSel	confTreePrev(ConfIn, NULL);
	cout << "OK." << endl;


INT	PrevMatDim = confTreePrev.init(0, NULL);
/*	for ( MRTreeIterator i = confTreePrev.firstInTree() ;
			!confTreePrev.isLastInTree(i) ; confTreePrev.nextInTree(i))
		PrevMatDim += confTreePrev.getConfigurationSAFNr(i).getSAFInc();
*/

	cout << "dimension of previous calculation: " << PrevMatDim << endl;

	cout << "building component mapping..." << flush;

	ticks.start();


INT	*Map = new INT[this->dim];
	for ( INT i=0 ; i<this->dim ; i++ )
		Map[i] = -1;

MRTreeIterator iPrev = confTreePrev.firstInTree();
ConfigurationSAFNr<MOType> confPrev = confTreePrev.getConfigurationSAFNr(iPrev);
INT	iSAF = 0;
const NExternalsDiag *confTree = extDiag;
	for ( MRTreeIterator i = confTree->firstInTree() ;
			!confTree->isLastInTree(i) ; confTree->nextInTree(i))
	{
	ConfigurationSAFNr<MOType> conf = confTree->getConfigurationSAFNr(i);
	INT	SAFs = conf.getSAFInc();

		if ( Pix p = ((extMOsSel *) iPrev.extMOs)->seek(conf) )
		{
		INT iSAFPrev = ((extMOsSel *) iPrev.extMOs)->
				getConfigurationSAFNr(p).getSAFNr();
				

			for ( INT k=0 ; k<SAFs ; k++ )
				Map[iSAF + k] = iSAFPrev++;

			if ( confTreePrev.isLastInTree(iPrev) )
				break;
			confTreePrev.nextInTree(iPrev);
			confPrev = confTreePrev.getConfigurationSAFNr(iPrev);
		}

		iSAF += SAFs;
	}

	ticks.stop();
	cout << ticks << " OK." << endl;



	cout << "projecting vectors to new space..." << flush;
	ticks.start();
	
Vector<CIVectorType>	ev(PrevMatDim);
DiskBuffer	evBuf(PrevMatDim*sizeof(CIVectorType),
	EigenvectorsFileName, DiskBuffer::noTempDir);

	for ( INT j=0 ; j<this->roots->getNumberOfRoots() ; j++ )
	{
		evBuf.get(j, ev.getP());
		memset(BasisVector, 0, this->dim*sizeof(CIVectorType));

		for ( INT k=0 ; k<this->dim ; k++ )
			if ( Map[k] != -1 )
				BasisVector[k] = ev[Map[k]];
		
//		for ( INT k=0 ; k<this->dim ; k++ )
//			cout << BasisVector[k] << endl;
		startVectors.put(j, BasisVector);
	}

	ticks.stop();
	cout << ticks << " OK." << endl;


	delete Map;
	delete BasisVector;

	start(&startVectors);
}

template <class CIMatrixType, class CIVectorType>
void	DavidsonCI<CIMatrixType, CIVectorType>::start(DiskBuffer *startVectors)
{
	Davidson<CIMatrixType, CIVectorType>::start(startVectors);

///////////////////////////////////////////////////////////////////////////////
//*
//*
	if ( method==ACPF || method==AQCC )
	{
	CIVectorType* p = new CIVectorType[this->dim];
		for ( INT i=0 ; i<this->basisDim ; i++ )
		{
			this->bBuf->get(i, p);
			calcExtPart(i+1, NULL);
//			calcExtPart(i+1, p);
		}
	delete p;
	}
//*
//*
///////////////////////////////////////////////////////////////////////////////
}


template <class CIMatrixType, class CIVectorType>
void	DavidsonCI<CIMatrixType, CIVectorType>::restart()
{
	Davidson<CIMatrixType, CIVectorType>::restart();
	
///////////////////////////////////////////////////////////////////////////////
//*
//*
	if ( method==ACPF || method==AQCC )
	{
	CIVectorType* p = new CIVectorType[this->dim];
		for ( INT i=0 ; i<this->basisDim ; i++ )
		{
			this->bBuf->get(i, p);
			if ( i<this->RefMatDim )
				calcExtPart(i+1, NULL);
			else
				calcExtPart(i+1, p);
		}
	delete p;
	}
//*
//*
///////////////////////////////////////////////////////////////////////////////
}




template <class CIMatrixType, class CIVectorType>
INT	DavidsonCI<CIMatrixType, CIVectorType>::doIteration(
	typename DavidsonCI<CIMatrixType, CIVectorType>::ConvergenceCriterion crit,
	double precision,
	INT &nconv, INT iter)
{
	if ( method==CI )
		return Davidson<CIMatrixType, CIVectorType>::doIteration(crit, precision, nconv, iter);
	else
		return doIterationSzalay(crit, precision, nconv, iter);
}


template <class CIMatrixType, class CIVectorType>
INT	DavidsonCI<CIMatrixType, CIVectorType>::doIterationSzalay(
	typename DavidsonCI<CIMatrixType, CIVectorType>::ConvergenceCriterion crit,
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
Vector<CIVectorType>	*y = new Vector<CIVectorType>(this->dim);


INT	bInc = 0;
DiskBuffer	*yBuf = new DiskBuffer(this->dim*sizeof(CIVectorType), 
	"Davidson_y.dat", DiskBuffer::truncateOnOpen);

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




///////////////////////////////////////////////////////////////////////////////
//*
//*	size consistency correcture
//*	see: P. Szalay
//*




	//	estimated correlation Energy
CIVectorType	*ECorr = new CIVectorType[this->roots->getNumberOfRoots()];

	//	weighting factor for ACPF, AQCC, ...
CIVectorType	g = 0;

INT	nElectrons = extDiag->getConfigurationSAFNr(0).getNumberOfElectrons();

//	cout << "electrons " << nElectrons << endl;

	switch ( method ) {
	case CI:
		g = 0;
		break;

	case ACPF:
		g = 1.0*(nElectrons-2)/nElectrons;
		break;
		
	case AQCC:
		g = 1.0*(nElectrons-2)*(nElectrons-3)/((nElectrons-1)*nElectrons);
		break;
	}



//*
//*
///////////////////////////////////////////////////////////////////////////////

	nconv = 0;

	/////////////////////////////////////////////////////////////////////
	//	optimized for minimum memory usage
	//


	for ( INT i=0 ; i<this->roots->getNumberOfRoots() ; i++ )
	{
		ECorr[i] = this->roots->getRefRoot(i)-this->lambda[this->roots->getRootNumber(i)];
//		cout << "!!" << this->roots->getRefRoot(i) << " " << this->lambda[this->roots->getRootNumber(i)] << endl;
	}


	for ( INT k=0 ; k<this->roots->getNumberOfRoots() ; k++ )
	{

///////////////////////////////////////////////////////////////////////////////
//*
//*
	CIMatrixType	*A = new CIMatrixType[this->basisDim*(this->basisDim+1)/2];

		memcpy(A, this->Aproj, this->basisDim*(this->basisDim+1)/2*sizeof(CIMatrixType));

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
		if ( ci2Ext )
			for ( INT i=0 ; i<this->basisDim ; i++ )
			{
				A[(i*i+3*i)/2] -= g*ECorr[k]*ci2Ext[i];
//					cout << "::" << ci2Ext[i] << "::" << endl;
			}
					
/*	cout << "projected matrix (shifted):" << endl;
	for ( INT i=0 ; i<basisDim ; i++ )
	{
		for ( INT j=0 ; j<basisDim ; j++ )
			if ( j<=i )
				cout << A[i*(i+1)/2 + j] << " ";
			else
				cout << "\t";
		cout << endl;
	}
	cout << endl;
*/

		hqrii1(
			&this->basisDim, &this->basisDim, A, this->lambda, this->alpha, &this->basisDim
		);

		delete A;

//	add some noise to start vectors to insure start of davidson algorithm
//	in pathologic cases
		if ( !iter )
		{
			if ( this->RefMatDim==1 )
			{
				for ( INT i=0 ; i<this->RefMatDim ; i++ )
						this->lambda[i] -= 0.01*(1.0*rand())/RAND_MAX;

				for ( INT i=0 ; i<this->RefMatDim ; i++ )
				{
					for ( INT j=0 ; j<this->RefMatDim ; j++ )
						this->alpha[j*this->RefMatDim + i] -= 0.01*(1.0*rand())/RAND_MAX;
				}
			}
		}



/*			cout << "projected eigenvalues:" << endl;
		for ( INT i=0 ; i<roots->getNumberOfRoots() ; i++ )
			cout << setw(22) << lambda[i] << endl;
		cout << endl;	
		cout << endl;

		cout << "projected eigenvectors (columnwise):" << endl;
		for ( INT i=0 ; i<basisDim ; i++ )
		{
			for ( INT j=0 ; j<basisDim ; j++ )
				cout << alpha[j*basisDim + i] << " ";
			cout << endl;
		}
		cout << endl;
*/		
//*
//*
///////////////////////////////////////////////////////////////////////////////


		if ( this->roots->getConvergenceStatus(k)==Roots::Converged )
		{
			nconv++;
			continue;
		}

		if ( crit == this->Energy )
			if ( fabs(this->lambda[this->roots->getRootNumber(k)] -
					this->oldLambda[this->roots->getRootNumber(k)])<precision )
			{
				this->roots->setConvergenceStatus(k, Roots::Converged);
				nconv++;
				continue;
			}


		x = new Vector<CIVectorType>(this->dim);
	
	CIMatrixType	lambdaK = this->lambda[this->roots->getRootNumber(k)];
	
		
	CIMatrixType	offset = g*ECorr[k];
//		cout << "offset=" << offset << endl;
		
		y->clear();
		
		for ( INT i=0 ; i<this->basisDim ; i++ )
		{
		CIMatrixType	alphaKI = this->alpha[this->roots->getRootNumber(k)*this->basisDim + i];

			{
			// calculate alpha_k*A*b_i
				this->AbBuf->get(i, x->getP());
		
			CIMatrixType	*pADiag = this->diagonal;
			CIVectorType	*pAb = x->getP();
			CIVectorType	*pf = y->getP();

//		cout << "lambdaK=" << lambdaK << endl;
//		cout << "*pADiag=" << *pADiag << endl;
//		cout << "alphaKI=" << alphaKI << endl;

				for ( INT j=0 ; j<internalCSFs ; j++ )
					if ( inactiveHole[j] )
						*pf++ += alphaKI/(lambdaK - *pADiag++ - offset) * *pAb++;
					else
						*pf++ += alphaKI/(lambdaK - *pADiag++) * *pAb++;

				for ( INT j=internalCSFs ; j<this->dim ; j++ )
					*pf++ += alphaKI/(lambdaK - *pADiag++ - offset) * *pAb++;
			}
			{
			// calculate alpha_k*lambda_k*b_i
				this->bBuf->get(i, x->getP());

			CIMatrixType	*pADiag = this->diagonal;
			CIVectorType	*pb = x->getP();
			CIVectorType	*pf = y->getP();

				for ( INT j=0 ; j<internalCSFs ; j++ )
					if ( inactiveHole[j] )
						*pf++ -= alphaKI/(lambdaK - *pADiag++ - offset) * lambdaK * *pb++;
					else
						*pf++ -= alphaKI/(lambdaK - *pADiag++) * lambdaK * *pb++;

				for ( INT j=internalCSFs ; j<this->dim ; j++ )
					*pf++ -= alphaKI/(lambdaK - *pADiag++ - offset) * lambdaK * *pb++;
			}
		}
		delete x;
		if ( crit == this->CorrectionVector )
			if ( fabs(y->getNorm2())<precision )
			{
				this->roots->setConvergenceStatus(k, Roots::Converged);
				nconv++;
				continue;
			}

		this->oldLambda[this->roots->getRootNumber(k)] = this->lambda[this->roots->getRootNumber(k)];

	
		y->normalize2();
		yBuf->put(k, y->getP());
	}
		
//	orthogonalize new vectors
	orthogonalize(this->bBuf, yBuf);
		
	for ( INT k=0 ; k<this->roots->getNumberOfRoots() ; k++ )
	{
		yBuf->get(k, y->getP());

//		printOrtho(bBuf);

//	check if new vector is worth to be appended to the basis vectors
//	cout << *y << endl;
	double	t = y->getNorm2();
//	cout << "Norm=" << t << endl;
		if ( t>this->SelThreshold )
		{
			*y /= t;
			this->bBuf->put(this->basisDim + bInc++, y->getP());

///////////////////////////////////////////////////////////////////////////////
//*
//*
			calcExtPart(this->basisDim + bInc, y->getP());
//*
//*
///////////////////////////////////////////////////////////////////////////////
		}
	}
	delete y;
	delete yBuf;

	delete ECorr;


	if ( !bInc )
		return nconv<this->roots->getNumberOfRoots();
	
	
//	calculate A*b
	this->davidsonMatrix->mult(this->bBuf, this->AbBuf, this->basisDim, this->basisDim+bInc-1);

	this->basisDim += bInc;

//	cout << bInc << " vectors appended to basis." << endl;
//	cout << "The basis now consists of " << basisDim << " vectors." << endl;



	allocateProjected(this->basisDim);

	calcProjected(this->basisDim-bInc);

	return 0;
}	







template <class CIMatrixType, class CIVectorType>
void	DavidsonCI<CIMatrixType, CIVectorType>::calcExtPart(INT BasisDim, CIVectorType *p)
{

	ci2Ext = (CIVectorType *) realloc(ci2Ext, BasisDim*sizeof(CIVectorType));

	ci2Ext[BasisDim-1] = 0;

	if ( p )
	{
	CIVectorType	h = 0;
	
		for ( INT i=0 ; i<internalCSFs ; i++ )
			if ( inactiveHole[i] )
				h += *p * *p++;
		ci2Ext[BasisDim-1] += h;
//		cout << "h0=" << h << endl;
		h = 0;
		for ( INT i=internalCSFs ; i<this->dim ; i++ )
			h += *p * *p++;
		ci2Ext[BasisDim-1] += h;
//		cout << "h1=" << h << endl;
	}
//	cout << "BasisDim-1=" << BasisDim-1 << endl;
}


template <class CIMatrixType, class CIVectorType>
INT		DavidsonCI<CIMatrixType, CIVectorType>::checkInactiveHole(
	Configuration<MOType> &conf,
	const MRMOs *mrmos)
{
	for ( INT i=0 ; i<mrmos->getNumberOfInactiveMOs() ; i++ )
	{
	INT	flag = 0;
/*		for ( INT j=0 ; j<conf.getNumberOfOpenShells() ; j++ )
			if ( mrmos->getInactiveMO(i) == conf.getOpenShell(j) )
				continue;
*/
		for ( INT j=0 ; j<conf.getNumberOfClosedShells() ; j++ )
			if ( mrmos->getInactiveMO(i) == conf.getClosedShell(j) )
			{
				flag = 1;
				continue;
			}

		if ( !flag )
			return 1;		
	}
	return 0;
}

template <class CIMatrixType, class CIVectorType>
void    DavidsonCI<CIMatrixType, CIVectorType>::iterate(
        INT maxIters,
        ConvergenceCriterion crit,
        double precision)
{
INT     nconv;
        for ( INT k=0 ; k<roots->getNumberOfRoots() ; k++ )
                roots->setConvergenceStatus(k, Roots::NotConverged);

        cout << setw(5) << 0;   
        for ( INT j=0 ; j<roots->getNumberOfRoots() ; j++ )
                cout << setw(22) << "---" <<
                        (roots->getConvergenceStatus(j)==Roots::Converged) ? " (*) " : "     ";
        cout << endl;   

        for ( INT i=0 ; i<roots->getNumberOfRoots() ; i++ )
                lambda[i] = roots->getRefRoot(i);
        
INT     status = 0;
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


//              cout << "iteration " << i+1 << " completed" << endl;
//              cout << nconv << " out of ";
//              cout << roots->getNumberOfRoots() << " roots converged" << endl;
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

//      if ( status==1 )
        {
        //      calculate eigenvectors
Vector<CIVectorType>       *ev = new Vector<CIVectorType>[roots->getNumberOfRoots()];
for ( INT i=0; i<roots->getNumberOfRoots(); i++)
{       
        ev[i].setDim(dim);
}
Vector<CIVectorType>    *bv = new Vector<CIVectorType>(dim);
DiskBuffer      *evBuf = new DiskBuffer(dim*sizeof(CIVectorType), "Eigenvectors.dat",
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


template class DavidsonCI<float, float>;
//template class DavidsonCI<float, double>;
//template class DavidsonCI<double, float>;
template class DavidsonCI<double, double>;
