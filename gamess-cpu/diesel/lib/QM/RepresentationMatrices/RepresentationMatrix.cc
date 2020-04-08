//***********************************************************************
//
//	Name:			RepresentationMatrix.cc
//
//	Description:			
//
//	Author:			Michael Hanrath
//
//	Version:		0.0
//
//	Date:			27.09.1996
//
//
//
//
//
//***********************************************************************

#include "RepresentationMatrix.h"

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <float.h>

using namespace std;

template <class MatrixType>
RepresentationMatrix<MatrixType>::RepresentationMatrix()
{
//	printf("###A %x\n", this);
//	fflush(stdout);
	elements = NULL;
	mpos = NULL;
	nSparse = nElements = 0;
	nMats = 1;
//	cout << "RepresentationMatrix<MatrixType>::RepresentationMatrix()" << endl;
//	cout << "MultTab[4]= " <<MultTab[4] << endl;
}


template <class MatrixType>
RepresentationMatrix<MatrixType>::RepresentationMatrix(
	MatrixType *buf, INT r, INT c, INT simplify)
{
//	printf("###A %x\n", this);
//	fflush(stdout);
	nMats = 1;
	round(buf, r*c);

	matDim.rows = r;
	matDim.cols = c;
INT	zero;
MatrixType	*p = buf;

/*	cout << "cols: " << matDim.cols << ", rows: " << matDim.rows << endl;

	for ( INT i=0 ; i<r ; i++ )
	{	for ( INT j=0 ; j<c ; j++ )
			cout << *p++ << " ";
		cout << endl;
	}
	p = buf;
*/

INT	isUnity = (matDim.rows==matDim.cols);
INT	isDiag = isUnity;
INT isTriDiag = isUnity;
INT	zeros = 0;

	elements = NULL;
	mpos = NULL;
	
	if (  matDim.rows*matDim.cols>=4 && simplify )
	{
		for ( INT i=0 ; i<matDim.rows ; i++ )
			for ( INT j=0 ; j<matDim.cols ; j++ , p++ )
			{	zero = (*p==0);
				if ( zero )
					zeros++;
				if ( (i!=j) && !zero )
					isUnity = isDiag = 0;
				if ( i==j && *p!=1.0 )
					isUnity = 0;
				if ( abs(i-j)>1 && !zero )
					isTriDiag = 0;
			}

//		isUnity = 0;
//		isDiag = 0;
//		isTriDiag = 0;
//		zeros = 0;

		if ( isUnity )
		{	elements = NULL;
			dense = unity;
			return;
		}

		if ( isDiag )
		{	elements = new MatrixType[matDim.cols];
			p = buf;
		MatrixType *pe = elements;
			for ( INT i=0 ; i<matDim.cols ; i++ , p+=matDim.cols+1)
				*pe++ = *p;
			dense = diagonal;
			nElements = matDim.cols;
			return;
		}

		if ( isTriDiag )
		{	elements = new MatrixType[3*matDim.cols-2];
			p = buf;
		MatrixType *pe = elements;
			*pe++ = *p++;
			*pe++ = *p++;
			for ( INT i=1 ; i<matDim.cols-1 ; i++ )
			{
				p+=matDim.cols-2;
				*pe++ = *p++;
				*pe++ = *p++;
				*pe++ = *p++;
			}
			p+=matDim.cols-2;
			*pe++ = *p++;
			*pe++ = *p++;
			dense = tridiagonal;
			nElements = 3*(matDim.cols-1) + 1;
			return;
		}

		if ( 8*zeros>matDim.rows*matDim.cols )
		{	nSparse = matDim.rows*matDim.cols - zeros;
			elements = new MatrixType[nSparse];
			mpos = new TMPos[nSparse];
			p = buf;
		MatrixType *pe = elements;

		INT	k = 0;
			for ( INT i=0 ; i<matDim.rows ; i++ )
				for ( INT j=0 ; j<matDim.cols ; j++ , p++ )
					if ( *p!=0 )
					{	mpos[k].r = i;
						mpos[k++].c = j;
						*pe++ = *p;
					}	

			dense = sparse;
			nElements = nSparse;
			return;
		}
	}

	elements = new MatrixType[matDim.rows*matDim.cols];
	memcpy(elements, buf, matDim.rows*matDim.cols*sizeof(MatrixType));
	dense = full;
	nElements = matDim.rows*matDim.cols;
	return;
}


template <class MatrixType>
RepresentationMatrix<MatrixType>::RepresentationMatrix(
	INT n, 
	MatrixType *buf, 
	INT r, INT c,
	MatrixDense _dense
	)
{
//	printf("###A %x\n", this);
//	fflush(stdout);
	nMats = n;
	switch ( _dense )
	{
	case full:
		nElements = n*r*c;
		break;

	case lowerTriangle:
		nElements = n*r*(r+1)/2;
		break;

	default:
		cerr << "error: wrong matrix type in RepresentationMatrix<MatrixType>::RepresentationMatrix" << endl;
	}
	round(buf, nElements);

	matDim.rows = r;
	matDim.cols = c;

	mpos = NULL;
	
	elements = new MatrixType[nElements];
	if ( elements == NULL )
        {
           cerr << "Not enough memory available for new MatrixType" << endl;
        }
	memcpy(elements, buf, nElements*sizeof(MatrixType));
	dense = _dense;
	return;
}


template <class MatrixType>
RepresentationMatrix<MatrixType>::RepresentationMatrix(const RepresentationMatrix &r)
{
//	printf("###A %x\n", this);
//	fflush(stdout);
	dense = r.getDense();
	matDim.rows = r.getNumberOfRows();
	matDim.cols = r.getNumberOfColumns();
	nSparse = r.getNSparse();
	nElements = r.getNElements();
	nMats = r.getNMats();
	mpos = NULL;
	elements = NULL;
	if ( nElements>0 )
	{
		elements = new MatrixType[nElements];
		memcpy(elements, r.getP(), nElements*sizeof(MatrixType));
		if ( dense==sparse  || dense==sameSparse )
		{
			mpos = new TMPos[nSparse];
			memcpy(mpos, r.getMPosP(), nSparse*sizeof(TMPos));
		}
	}
}


template <class MatrixType>
RepresentationMatrix<MatrixType>::~RepresentationMatrix()
{
//	printf("###D %x\n", this);
//	fflush(stdout);
	if ( elements )
		delete[] elements;
	if ( mpos )
		delete[] mpos;
}


template <class MatrixType>
void	RepresentationMatrix<MatrixType>::round(MatrixType *p, INT n)
{
double	eps = 1e6*DBL_EPSILON;
	for ( INT i=0 ; i<n ; i++ , p++ )
	{	if ( fabs(*p)<eps )
			*p = 0;
		if ( fabs(*p - 1.0)<eps )
			*p = 1.0;
	}
}


template <class MatrixType>
INT	RepresentationMatrix<MatrixType>::getOccupiedMemory() const
{
INT	l = sizeof(RepresentationMatrix);

	switch ( dense ) {
	case unity:
		break;

	case diagonal:
		l += matDim.cols*sizeof(MatrixType);
		break;
		
	case tridiagonal:
		l += (3*matDim.cols-2)*sizeof(MatrixType);
		break;
		
	case sparse:
	case sameSparse:
		l += nSparse*(sizeof(MatrixType)+2*sizeof(INT));
		break;
		
	case full:
		l += matDim.cols*matDim.rows*sizeof(MatrixType);
		break;

	case lowerTriangle:
		l += matDim.cols*(matDim.cols+1)/2*sizeof(MatrixType);
		break;
	}
	return l;
}

//-----------------------------------------------------------------------
//-----------------------------------------------------------------------
//-----------------------------------------------------------------------

template <class MatrixType>
ostream &	operator<<(ostream &s, const RepresentationMatrix<MatrixType> & repMat)
{
const MatrixType	*p = repMat.getP();

	switch ( repMat.getDense() ) {
	case RepresentationMatrix<MatrixType>::unity:
		s << "unity" << endl;
		for ( INT i=0 ; i<repMat.getNumberOfRows() ; i++ )
		{	s << "(";
			for ( INT j=0 ; j<repMat.getNumberOfColumns() ; j++ )
			{	if ( i==j )
					s << 1;
				else
					s << 0;
				if ( j<repMat.getNumberOfColumns()-1 )
					s << " ";
			}
			s << ")" << endl;
		}
		break;

	case RepresentationMatrix<MatrixType>::diagonal:
		s << "diagonal" << endl;
		for ( INT i=0 ; i<repMat.getNumberOfRows() ; i++ )
		{	s << "(";
			for ( INT j=0 ; j<repMat.getNumberOfColumns() ; j++ )
			{	if ( i==j )
					s << *p++;
				else
					s << 0;
				if ( j<repMat.getNumberOfColumns()-1 )
					s << " ";
			}
			s << ")" << endl;
		}
		break;
		
	case RepresentationMatrix<MatrixType>::tridiagonal:
		s << "tridiagonal" << endl;
		for ( INT i=0 ; i<repMat.getNumberOfRows() ; i++ )
		{	s << "(";
			for ( INT j=0 ; j<repMat.getNumberOfColumns() ; j++ )
			{	if ( abs(i-j)<=1 )
					s << *p++;
				else
					s << 0;
				if ( j<repMat.getNumberOfColumns()-1 )
					s << " ";
			}
			s << ")" << endl;
		}
		break;
		
	case RepresentationMatrix<MatrixType>::sparse:
		s << "sparse" << endl;
		for ( INT i=0 ; i<repMat.getNumberOfRows() ; i++ )
		{	s << "(";
			for ( INT j=0 ; j<repMat.getNumberOfColumns() ; j++ )
			{	
			INT found = 0;
				for ( INT k=0 ; k<repMat.getNSparse() ; k++ )
					if ( repMat.getMPosP()[k].r==i && repMat.getMPosP()[k].c==j )
					{	found = 1;
						break;
					}
				if ( found )
					s << *p++;
				else
					s << "0";
				if ( j<repMat.getNumberOfColumns()-1 )
					s << " ";
			}
			s << ")" << endl;
		}

		break;
		
	case RepresentationMatrix<MatrixType>::lowerTriangle:
		s << "lowerTriangle" << endl;
		cout <<repMat.nMats << " " << repMat.getNumberOfRows()<< " " << repMat.getNumberOfRows()<< endl;
		for ( INT k=0 ; k<repMat.nMats ; k++ )
		{
			for ( INT i=0 ; i<repMat.getNumberOfRows() ; i++ )
			{	s << "(";
				for ( INT j=0 ; j<repMat.getNumberOfRows() ; j++ )
				{
					if ( j<=i )
						s << p[(i+1)*i/2*repMat.nMats +
							j*repMat.nMats + k] << " ";
					else
						s << "\t" << " ";
				}
				s << ")" << endl;
			}
			cout << endl;
		}
		break;
		
	case RepresentationMatrix<MatrixType>::full:
		s << "full" << endl;
		for ( INT k=0 ; k<repMat.nMats ; k++ )
		{
			for ( INT i=0 ; i<repMat.getNumberOfRows() ; i++ )
			{	s << "(";
				for ( INT j=repMat.getNumberOfColumns()-1 ; j>=0 ; j-- )
				{	
					s << p[(i+1)*i/2*repMat.nMats +
						j*repMat.nMats + k] << " ";
					if ( j )
						s << " ";
				}
				s << ")" << endl;
			}
			cout << endl;
		}
		break;

	case RepresentationMatrix<MatrixType>::sameSparse:
		s << "sameSparse" << endl;
		break;
	}
	

	return s;
}


template class RepresentationMatrix<float>;
template class RepresentationMatrix<double>;

template ostream & operator << (ostream & s, const RepresentationMatrix<double> &);	
template ostream & operator << (ostream & s, const RepresentationMatrix<float> &);	
