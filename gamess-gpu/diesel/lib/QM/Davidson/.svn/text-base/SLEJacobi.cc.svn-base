//***********************************************************************
//
//	Name:			SLEJacobi.cc
//
//	Description:	solution of large linear equation system
//					by Jacobi's method 
//					(total step method: 
//						inverse of A is approximated
//						by inverse of diagonal)
//
//	Author:			Michael Hanrath
//
//	Version:		0.0
//
//	Date:			16.08.1998
//
//
//
//
//***********************************************************************


#include "SLEJacobi.h"
#include "SLEMatrix.h"
#include "../../Math/MatrixVector/BufferedVector.h"

#include "../../Math/MatrixVector/Vector.h"

#include <math.h>

using namespace std;

template <class MatrixType, class VectorType>
SLEJacobi<MatrixType, VectorType>::SLEJacobi(
	const SLEMatrix<MatrixType, VectorType> *_sleMatrix) :
	x(*(new BufferedVector<VectorType>(_sleMatrix->getTotalDim()))),
	b(*(new BufferedVector<VectorType>(_sleMatrix->getTotalDim())))
{
	sleMatrix = _sleMatrix;
}


template <class MatrixType, class VectorType>
SLEJacobi<MatrixType, VectorType>::~SLEJacobi()
{
	delete &x;
	delete &b;
}


template <class MatrixType, class VectorType>
const	BufferedVector<VectorType>	&SLEJacobi<MatrixType, VectorType>::getX() const
{
	return	x;
}

template <class MatrixType, class VectorType>
const	BufferedVector<VectorType>	&SLEJacobi<MatrixType, VectorType>::getB() const
{
	return	b;
}


template <class MatrixType, class VectorType>
void	SLEJacobi<MatrixType, VectorType>::iterate()
{
//	Problem:
//
//		solve A*x=b for x, A: VERY large matrix
//
//
//	Method (Jacobi, total step method):
//
//                  1
// 1.       x  = -------  b 
//           i      A      i
//                   ii
//
// 2.       r  = b - A*x
//
//                  1
// 3.       d  = -------  r 
//           i      A      i
//                   ii
//
// 4.       x  = x + d
//
// 5.       if  ||d|| > thresh goto 2.
//
//

BufferedVector<VectorType>	y(sleMatrix->getTotalDim());

/*

	for ( INT i=0 ; i<sleMatrix->getTotalDim() ; i++ )
	{
		cout << "=:=:=:= " << i << "=:=:=:=" << endl;
		vx.clear();
		vx[i] = 1;
		xBuf->put(0, vx.getP());
		vy.clear();
		yBuf->put(0, vy.getP());
		sleMatrix->mult(xBuf, yBuf);
		yBuf->get(0, vy.getP());
		cout << vy << endl;
	}



	return;


*/

	sleMatrix->calcResidual(x);

VectorType	normb = x.getInfNorm();
	b = x;

	sleMatrix->multInvDiag(x);
INT	iter = 0;
	while ( 1 )
	{
		cout << "iteration #" << ++iter << endl;
		sleMatrix->mult(x, y);
		y -= b;
		sleMatrix->multInvDiag(y);
		
		x -= y;

VectorType	norm = y.getInfNorm();
		cout << "infNorm=" << norm << endl;

		if ( fabs(norm/normb)<epsilon )
			break;

	}
/*	cout << "Probe: " << endl;
	
	sleMatrix->mult(xBuf, yBuf);
	yBuf->get(0, vy.getP());
	cout << b << endl;
	cout << vy << endl;
*/
}






template class SLEJacobi<float, float>;
template class SLEJacobi<double, float>;
template class SLEJacobi<float, double>;
template class SLEJacobi<double, double>;

template <class MatrixType, class VectorType>
const VectorType SLEJacobi<MatrixType, VectorType>::epsilon = 1e-4;
