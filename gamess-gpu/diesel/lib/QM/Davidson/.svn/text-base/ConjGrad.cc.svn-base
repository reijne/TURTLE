//***********************************************************************
//
//	Name:			ConjGrad.cc
//
//	Description:	solution of large linear equation system
//					by method of conjugated gradient
//
//	Author:			Michael Hanrath
//
//	Version:		0.0
//
//	Date:			30.08.1998
//
//	Ref:			
//
//
//***********************************************************************


#include "ConjGrad.h"
#include "SLEMatrix.h"

#include "../../Math/MatrixVector/BufferedVector.h"

#include <math.h>

#include <stdlib.h>

using namespace std;

void	lst()
{
	cout << "==============================" << endl;
	system("ls -l Vector*");
	cout << "==============================" << endl;
	cout << endl;
	cout << endl;
	cout << endl;
}

template <class MatrixType, class VectorType>
ConjGrad<MatrixType, VectorType>::ConjGrad(
	const SLEMatrix<MatrixType, VectorType> *_sleMatrix) :
	x(*(new BufferedVector<VectorType>(_sleMatrix->getTotalDim()))),
	b(*(new BufferedVector<VectorType>(_sleMatrix->getTotalDim())))
{
	sleMatrix = _sleMatrix;
}


template <class MatrixType, class VectorType>
ConjGrad<MatrixType, VectorType>::~ConjGrad()
{
	delete &x;
	delete &b;
}


template <class MatrixType, class VectorType>
const	BufferedVector<VectorType>	&ConjGrad<MatrixType, VectorType>::getX() const
{
	return	x;
}


template <class MatrixType, class VectorType>
const	BufferedVector<VectorType>	&ConjGrad<MatrixType, VectorType>::getB() const
{
	return	b;
}


template <class MatrixType, class VectorType>
void	ConjGrad<MatrixType, VectorType>::iterate()
{
INT	totalDim = sleMatrix->getTotalDim();


BufferedVector<VectorType>	p(totalDim);
BufferedVector<VectorType>	q(totalDim);
BufferedVector<VectorType>	r(totalDim);
BufferedVector<VectorType>	z(totalDim);


	sleMatrix->calcResidual(x);
	x = -x;
	b = x;

	cout << "bnorm2= " << b.getNorm2Sqr() << endl;

	{
	BufferedVector<VectorType>	y(totalDim);
		sleMatrix->multInvDiag(x);
		sleMatrix->mult(x, y);
		r = b-y;
	}

VectorType	normb = x.getInfNorm();
	if ( fabs(normb) < epsilon )
	{
		cout << "residual norm to small" << endl;
	}
	else
	{

	VectorType	rho = 0;
	VectorType	rho1 = 0;
	VectorType	alpha = 0;
	VectorType	beta = 0;



	INT	iter = 0;
		while ( iter<50 )
		{
			cout << "iteration #" << ++iter << endl;

	//		z = D*r;
			z = r;
			sleMatrix->multInvDiag(z);

			rho = r*z;

			if ( iter==1 )
				p = z;
			else
			{
				beta = rho/rho1;
				p = z + beta*p;
			}

	//		q = A*p;
			sleMatrix->mult(p, q);

			alpha = rho/(p*q);

			x += alpha*p;
			r -= alpha*q;
			

		VectorType	resid = r.getInfNorm();
			cout << "infNorm=" << resid << endl;

			if ( resid/normb<epsilon )
				break;

			rho1 = rho;
		}
	}
	
/*	cout << "x = " << x << endl;
	cout << "Probe: " << endl;
	
	sleMatrix->mult(x, y);
	cout << b << endl;
	cout << y << endl;
*/
}



/*
template <class MatrixType, class VectorType>
void	ConjGrad<MatrixType, VectorType>::iterate()
{
INT	totalDim = sleMatrix->getTotalDim();
BufferedVector<VectorType>	&y = *(new BufferedVector<VectorType>(totalDim));


BufferedVector<VectorType>	p(totalDim);
BufferedVector<VectorType>	q(totalDim);
BufferedVector<VectorType>	r(totalDim);
BufferedVector<VectorType>	z(totalDim);


	sleMatrix->calcResidual(x);
	x = -x;
	b = x;

	cout << "bnorm2= " << b.getNorm2Sqr() << endl;

	sleMatrix->multInvDiag(x);
	sleMatrix->mult(x, y);
	r = b-y;

VectorType	normb = x.getInfNorm();
	if ( fabs(normb) < epsilon )
	{
		cout << "residual norm to small" << endl;
	}
	else
	{

	VectorType	rho = 0;
	VectorType	rho1 = 0;
	VectorType	alpha = 0;
	VectorType	beta = 0;



	INT	iter = 0;
		while ( iter<50 )
		{
			cout << "iteration #" << ++iter << endl;


	//		z = D*r;
			z = r;
			sleMatrix->multInvDiag(z);

			rho = r*z;

			if ( iter==1 )
				p = z;
			else
			{
				beta = rho/rho1;
				p = z + beta*p;
			}

	//		q = A*p;
			sleMatrix->mult(p, q);

			alpha = rho/(p*q);

			x += alpha*p;
			r -= alpha*q;
			

		VectorType	resid = r.getInfNorm();
			cout << "infNorm=" << resid << endl;

			if ( resid/normb<epsilon )
				break;

			rho1 = rho;
		}
	}
	
//	cout << "x = " << x << endl;
//	cout << "Probe: " << endl;
	
//	sleMatrix->mult(x, y);
//	cout << b << endl;
//	cout << y << endl;

	delete &y;
}


*/



template <class MatrixType, class VectorType>
void	ConjGrad<MatrixType, VectorType>::iterateS()
{
/*
INT	totalDim = sleMatrix->getTotalDim();

BufferedVector<VectorType>	&y = *(new BufferedVector<VectorType>(totalDim));



BufferedVector<VectorType>	p(totalDim);
BufferedVector<VectorType>	phat(totalDim);
BufferedVector<VectorType>	q(totalDim);
BufferedVector<VectorType>	qhat(totalDim);
BufferedVector<VectorType>	u(totalDim);
BufferedVector<VectorType>	uhat(totalDim);
BufferedVector<VectorType>	vhat(totalDim);
BufferedVector<VectorType>	r(totalDim);
BufferedVector<VectorType>	rtilde(totalDim);
BufferedVector<VectorType>	h(totalDim);



	sleMatrix->calcResidual(xBuf);
	xBuf->get(0, x.getP());
	p = -x;
	xBuf->put(0, p.getP());
//	cout << "b= " << x << endl;

VectorType	normb = x.getInfNorm();
	if ( fabs(normb) < epsilon )
	{
		cout << "residual norm to small" << endl;
		exit(1);
	}

	sleMatrix->multInvDiag(xBuf);
	xBuf->get(0, x.getP());
//	cout << "b(inv)= " << x << endl;

	sleMatrix->mult(xBuf, yBuf);
	sleMatrix->calcResidual(yBuf);
	yBuf->get(0, r.getP());
	
	rtilde = r;
	
VectorType	rho1 = 0;
VectorType	rho2 = 0;
VectorType	alpha = 0;
VectorType	beta = 0;



INT	iter = 0;
	while ( 1 )
	{
		cout << "iteration #" << ++iter << endl;
		
		rho1 = rtilde*r;
		if ( fabs(rho1) == 0 )
		{
			cout << "residual norm to small" << endl;
			exit(1);
		}
		
		if ( iter==1 )
		{
			u = r;
			p = u;
		}
		else
		{
			beta = rho1/rho2;
			u = r + beta*q;
			p = u + beta*(q + beta*p);
		}
		
//		phat = D*p;
		xBuf->put(0, p.getP());
		sleMatrix->multInvDiag(xBuf);
		xBuf->get(0, phat.getP());
		
		
//		vhat = A*phat;
		xBuf->put(0, phat.getP());
		sleMatrix->mult(xBuf, yBuf);
		yBuf->get(0, vhat.getP());
		
		alpha = rho1/(rtilde*vhat);
		
		q = u - alpha*vhat;
		


//		uhat = D*(u+q);
		h = u + q;
		xBuf->put(0, h.getP());
		sleMatrix->multInvDiag(xBuf);
		xBuf->get(0, uhat.getP());
		
		x += alpha*uhat;
		
//		qhat = A*uhat;
		xBuf->put(0, uhat.getP());
		sleMatrix->mult(xBuf, yBuf);
		yBuf->get(0, qhat.getP());
		
		r -= alpha*qhat;
		
		rho2 = rho1;
				
	VectorType	resid = r.getInfNorm();
		cout << "infNorm=" << resid << endl;
	
		if ( resid/normb<epsilon )
			break;

	}

	delete &y;
*/
}





template class ConjGrad<float, float>;
template class ConjGrad<double, float>;
template class ConjGrad<float, double>;
template class ConjGrad<double, double>;

template <class MatrixType, class VectorType>
const VectorType ConjGrad<MatrixType, VectorType>::epsilon=1e-4;
