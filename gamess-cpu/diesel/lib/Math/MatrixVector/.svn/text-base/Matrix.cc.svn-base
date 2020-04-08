//***********************************************************************
//
//	Name:			Matrix.cc
//
//	Description:	
//
//	Author:			Michael Hanrath
//
//	Version:		0.0
//
//	Date:			05. Nov 1998
//
//***********************************************************************




#include "Matrix.h"

#include "../FortranNumeric/FortranEigenProblems.h"


#include <string>
#include <stdlib.h>
#include <iomanip>
#include <streambuf>
#include <assert.h>

using namespace std;

template <>
Matrix<INT>::Matrix(INT _r, INT _c, INT a1 ...)
{
INT	i, j;
va_list	ap;
INT	*pl;

	r = _r;
	c = _c;
	p = new INT[r*c];
	if ( !p )
		Fehler(0);

	pl = p;
	*pl++ = a1;
	va_start(ap, a1);
	for ( i=0 ; i<r ; i++ )
		for ( j=0 ; j<c ; j++ )
			if ( i+j )
				*pl++ = va_arg(ap, INT);
	va_end(ap);
}

template <>
Matrix<float>::Matrix(INT _r, INT _c, INT a1 ...)
{
INT	i, j;
va_list	ap;
float	*pl;

	r = _r;
	c = _c;
	p = new float[r*c];
	if ( !p )
		Fehler(0);

	pl = p;
	*pl++ = a1;
	va_start(ap, a1);
	for ( i=0 ; i<r ; i++ )
		for ( j=0 ; j<c ; j++ )
			if ( i+j )
				*pl++ = va_arg(ap, INT);
	va_end(ap);
}

template <>
Matrix<float>::Matrix(INT _r, INT _c, double a1 ...)
{
INT	i, j;
va_list	ap;
float	*pl;

	r = _r;
	c = _c;
	p = new float[r*c];
	if ( !p )
		Fehler(0);

	pl = p;
	*pl++ = a1;
	va_start(ap, a1);
	for ( i=0 ; i<r ; i++ )
		for ( j=0 ; j<c ; j++ )
			if ( i+j )
				*pl++ = va_arg(ap, double);
	va_end(ap);
}

template <>
Matrix<double>::Matrix(INT _r, INT _c, double a1 ...)
{
INT	i, j;
va_list	ap;
double	*pl;

	r = _r;
	c = _c;
	p = new double[r*c];
	if ( !p )
		Fehler(0);

	pl = p;
	*pl++ = a1;
	va_start(ap, a1);
	for ( i=0 ; i<r ; i++ )
		for ( j=0 ; j<c ; j++ )
			if ( i+j )
				*pl++ = va_arg(ap, double);
	va_end(ap);
}


template <>
Matrix<double>::Matrix(INT _r, INT _c, INT a1 ...)
{
INT	i, j;
va_list	ap;
double	*pl;

	r = _r;
	c = _c;
	p = new double[r*c];
	if ( !p )
		Fehler(0);

	pl = p;
	*pl++ = a1;
	va_start(ap, a1);
	for ( i=0 ; i<r ; i++ )
		for ( j=0 ; j<c ; j++ )
			if ( i+j )
				*pl++ = va_arg(ap, INT);
	va_end(ap);
}

/*

template <>
Matrix<INT>::operator Matrix<double>()
{
INT	i;
double	*A;

	A = new double[r*c];
	for ( i=0 ; i<r*c ; i++ )
		*(A+i) = p[i];

	return Matrix<double>(r, c, A);
}


template <>
Matrix<float>::operator Matrix<double>()
{
INT	i;
Matrix<double>	*A;

	A = new Matrix<double>(r, c);
	for ( i=0 ; i<r*c ; i++ )
		(*A)[i] = p[i];

	return *A;
}


template <>
Matrix<double>::operator Matrix<Complex>()
{
INT	i;
Matrix<Complex>	*A;

	A = new Matrix<Complex>(r, c);
	for ( i=0 ; i<r*c ; i++ )
		(*A)[i] = p[i];

	return *A;
}
*/

/*

template <>
Matrix<Complex> Complex(Matrix<double> & vin)
{
Matrix<Complex>	*v;
INT	i;

	v = new Matrix<Complex>(vin.getRows(), vin.getCols());

	for ( i=0 ; i<vin.getRows() * vin.getCols() ; i++ )
		(*v)[i] = ::Complex(vin[i]);
	return *v;
}
*/

/*
Matrix<double> real(const Matrix<Complex> & vin)
{
Matrix<double>	*v;
INT	i;

	v = new Matrix<double>(vin.getRows(), vin.getCols());

	for ( i=0 ; i<vin.getRows() * vin.getCols() ; i++ )
		(*v)[i] = ::real(vin[i]);
	return *v;
}

Matrix<double> imag(const Matrix<Complex> & vin)
{
Matrix<double>	*v;
INT	i;

	v = new Matrix<double>(vin.getRows(), vin.getCols());

	for ( i=0 ; i<vin.getRows() * vin.getCols() ; i++ )
		(*v)[i] = ::imag(vin[i]);
	return *v;
}
*/

Matrix<INT> integer(const Matrix<double> & vin)
{
Matrix<INT>	v(vin.getRows(), vin.getCols());
INT	i;


	for ( i=0 ; i<vin.getRows() * vin.getCols() ; i++ )
		(v)[i] = (INT) (vin[i]+0.5);
	return v;
}





//-----------------------------------------------------------------------

Matrix<Complex> operator+(const Matrix<Complex> &v1, const Matrix<double> &v2)
{
INT	i;
	if ( (v1.getRows() != v2.getRows()) || (v1.getCols() != v2.getCols()) )
	{	v1.Fehler(1);
Matrix<Complex> 	v;
		return v;
	}
Matrix<Complex> 	v(v1.getRows(), v1.getCols());
	for ( i=0 ; i<v1.getRows()*v1.getCols() ; i++ )
		v[i] = v1[i] + v2[i];
	return v;
}

Matrix<Complex> operator-(const Matrix<Complex> &v1, const Matrix<double> &v2)
{
INT	i;
	if ( (v1.getRows() != v2.getRows()) || (v1.getCols() != v2.getCols()) )
	{	v1.Fehler(1);
Matrix<Complex> 	v;
		return v;
	}
Matrix<Complex>	v(v1.getRows(), v1.getCols());
	for ( i=0 ; i<v1.getRows()*v1.getCols() ; i++ )
		v[i] = v1[i] - v2[i];
	return v;
}


Matrix<Complex> operator*(const Matrix<Complex> &A, const Matrix<double> &B)
{
INT	p, q, r;
INT	i, j, k;
Complex	*pA, *pC;
double	*pB;
	if ( A.getCols() != B.getRows() )
	{	A.Fehler(1);
Matrix<Complex> 	C;
		return C;
	}
	p = A.getRows();
	q = A.getCols();
	r = B.getCols();

Matrix<Complex>	C(p, r);

	for ( i=0 ; i<p ; i++ )
		for ( j=0 ; j<r ; j++ )
		{   pA = &A[0] + i*q;
			pB = &B[0] + j;
			pC = &C[0] + i*r + j;
			*pC = 0;
			for ( k=0 ; k<q ; k++ )
			{	*pC += *pA * *pB;
				pA++;
				pB += r;
			}
		}

	return C;
}




Matrix<Complex> operator+(const Matrix<double> &v1,const  Matrix<Complex> &v2)
{	return	v2+v1;	};

Matrix<Complex> operator-(const Matrix<double> &v1, const Matrix<Complex> &v2)
{	return	-(v2-v1);	};

Matrix<Complex> operator*(const Matrix<double> &A, const Matrix<Complex> &B)
{
INT	p, q, r;
INT	i, j, k;
Complex	*pB, *pC;
double	*pA;
	if ( A.getCols() != B.getRows() )
	{	A.Fehler(1);
Matrix<Complex> 	C;
		return C;
	}
	p = A.getRows();
	q = A.getCols();
	r = B.getCols();

Matrix<Complex>	C(p, r);

	for ( i=0 ; i<p ; i++ )
		for ( j=0 ; j<r ; j++ )
		{   pA = &A[0] + i*q;
			pB = &B[0] + j;
			pC = &C[0] + i*r + j;
			*pC = 0;
			for ( k=0 ; k<q ; k++ )
			{	*pC += *pA * *pB;
				pA++;
				pB += r;
			}
		}

	return C;
}







//-----------------------------------------------------------------------


Vector<Complex> operator*(const Matrix<double> &A, const Vector<Complex> &v1)
{
INT	p, q;
INT	i, j;
double	*pA;
	if ( A.getCols() != v1.getDim() )
	{	A.Fehler(1);
Vector<Complex> 	v;
		return v;
	}
	p = A.getRows();
	q = A.getCols();

Vector<Complex>	v(p);

	pA = &A[0];

	for ( i=0 ; i<p ; i++ )
	{   v[i] = 0;
		for ( j=0 ; j<q ; j++ )
			v[i] += *pA++ * v1[j];
	}

	return v;
}

Vector<Complex> operator*(const Matrix<Complex> &A, const Vector<double> &v1)
{
INT	p, q;
INT	i, j;
Complex	*pA;
	if ( A.getCols() != v1.getDim() )
	{	A.Fehler(1);
Vector<Complex> 	v;
		return v;
	}
	p = A.getRows();
	q = A.getCols();

Vector<Complex>	v(p);

	pA = &A[0];

	for ( i=0 ; i<p ; i++ )
	{   v[i] = 0;
		for ( j=0 ; j<q ; j++ )
			v[i] += *pA++ * v1[j];
	}

	return v;
}




//-------------------------------------------------------------------
//*************************************************************************
//*************************************************************************
//*************************************************************************
//*************************************************************************




template <class T> Matrix<T>::Matrix(INT _r, INT _c)
{
	r = _r;
	c = _c;
	p = new T[r*c];
  	if ( !p )
		Fehler(0);
}

template <class T> Matrix<T>::Matrix(INT _r, INT _c, const T *_p)
{
	r = _r;
	c = _c;
	p = new T[r*c];
	memcpy(p, _p, r*c*sizeof(T));
}



template <class T> Matrix<T>::~Matrix()
{
	if ( r*c )
		delete p;
}


template <class T> Matrix<T>& Matrix<T>::operator=(const Matrix<T> & v1)
{
INT	i;

	if ( (v1.getRows() != r) || (v1.getCols() != c) )
	{	if ( r+c )
			delete	p;
		r = v1.getRows();
		c = v1.getCols();
		p = new T[r*c];
	}
	for ( i=0 ; i<r*c ; i++ )
		p[i] = v1[i];
	return *this;
}



template <class T> Matrix<T>::Matrix(const Matrix<T> & v)
{
	c = v.getCols();
	r = v.getRows();
	p = new T[r*c];
	memcpy(p, v.getP(), r*c*sizeof(T));
}








template <class T> Matrix<T> operator*(const T a, const Matrix<T> &A)
{
INT	p, q;
INT	i;
T	*pA, *pC;

	p = A.getRows();
	q = A.getCols();

Matrix<T>	C(p, q);

	pA = &A[0];
	pC = &C[0];

	for ( i=0 ; i<p*q ; i++ )
		*pC++ = *pA++ * a;

	return C;
}







template <class T> Matrix<T> operator+(const Matrix<T> &v1, const Matrix<T> &v2)
{
INT	i;
T 	*v;
	if ( (v1.getRows() != v2.getRows()) || (v1.getCols() != v2.getCols()) )
	{	v1.Fehler(1);
		return Matrix<T>();
	}
	v = new T[v1.getRows()*v1.getCols()];
	for ( i=0 ; i<v1.getRows()*v1.getCols() ; i++ )
		*(v+i) = v1[i] + v2[i];
	return Matrix<T>(v1.getRows(), v1.getCols(), v);
}

template <class T> Matrix<T> operator-(const Matrix<T> &v1, const Matrix<T> &v2)
{
INT	i;
	if ( (v1.getRows() != v2.getRows()) || (v1.getCols() != v2.getCols()) )
	{	v1.Fehler(1);
Matrix<T> 	v;
		return v;
	}
Matrix<T>	v(v1.getRows(), v1.getCols());
	for ( i=0 ; i<v1.getRows()*v1.getCols() ; i++ )
		v[i] = v1[i] - v2[i];
	return v;
}

template <class T> Matrix<T> operator*(const Matrix<T> &A, const Matrix<T> &B)
{
INT	p, q, r;
INT	i, j, k;
T 	*C;
T	*pA, *pB, *pC;

	if ( A.getCols() != B.getRows() )
	{	A.Fehler(1);
		return Matrix<T>();
	}
	p = A.getRows();
	q = A.getCols();
	r = B.getCols();

	C = new T[p*r];

	for ( i=0 ; i<p ; i++ )
		for ( j=0 ; j<r ; j++ )
		{   pA = &A[0] + i*q;
			pB = &B[0] + j;
			pC = C + i*r + j;
			*pC = 0;
			for ( k=0 ; k<q ; k++ )
			{	*pC += *pA * *pB;
				pA++;
				pB += r;
			}
		}

	return Matrix<T>(p, r, C);
}


template <class T> Vector<T> operator*(const Matrix<T> &A, const Vector<T> &v1)
{
INT	p, q;
INT	i, j;
T	*pA;
	if ( A.getCols() != v1.getDim() )
	{	A.Fehler(1);
Vector<T> 	v;
		return v;
	}
	p = A.getRows();
	q = A.getCols();

Vector<T>	v(p);

	pA = &A[0];

	for ( i=0 ; i<p ; i++ )
	{   v[i] = 0;
		for ( j=0 ; j<q ; j++ )
			v[i] += *pA++ * v1[j];
	}

	return v;
}

/*
template <class T> Vector<T> operator*(const Vector<T> &v1, const Matrix<T> &A)
{
	return A*v1;
}
*/


template <class T> Matrix<T> & Matrix<T>::operator+=(const Matrix<T> &v1)
{
INT	i;
	if ( (this->getRows() != v1.getRows()) || (this->getCols() != v1.getCols()) )
	{	v1.Fehler(1);
		return *this;
	}
	for ( i=0 ; i<v1.getRows() * v1.getCols() ; i++ )
		(*this)[i] += v1[i];
	return *this;
}

template <class T> Matrix<T> & Matrix<T>::operator-=(const Matrix<T> &v1)
{
INT	i;
	if ( (this->getRows() != v1.getRows()) || (this->getCols() != v1.getCols()) )
	{	v1.Fehler(1);
		return *this;
	}
	for ( i=0 ; i<v1.getRows() * v1.getCols() ; i++ )
		(*this)[i] -= v1[i];
	return *this;
}


template <class T> Matrix<T> & Matrix<T>::operator*=(const Matrix<T> &B)
{
INT	p, q, r;
INT	i, j, k;
Matrix<T> 	*A;
Matrix<T> *BB;
T		*pA, *pB, *pC;

	if ( this->getCols() != B.getRows() )
	{	this->Fehler(1);
		return *this;
	}
	p = this->getRows();
	q = this->getCols();
	r = B.getCols();

	A = new Matrix<T>(p, r);
	*A = *this;
	if ( &(*this)[0] == &B[0] )
	{	BB = new Matrix<T>(q, r);
		*BB = B;
	}
	else
		BB = (Matrix<T> *) &B;

	for ( i=0 ; i<p ; i++ )
		for ( j=0 ; j<r ; j++ )
		{   pA = &(*A)[0] + i*q;
			pB = &(*BB)[0] + j;
			pC = &(*this)[0] + i*r + j;
			*pC = 0;
			for ( k=0 ; k<q ; k++ )
			{	*pC += *pA * *pB;
				pA++;
				pB += r;
			}
		}
	delete A;
	if ( &(*this)[0] == &B[0] )
    	delete BB;
	return *this;
}



template <class T> Matrix<T> & Matrix<T>::operator*=(const T & mul)
{
T	*pA = &(*this)[0];
	for ( INT i=0 ; i<r*c ; i++ )
		*pA++ *= mul;
		
	return *this;
}



template <class T> void Matrix<T>::setFill(const T & value)
{	
T	*pA = &(*this)[0];

	for ( INT i=0 ; i<r*c ; i++ )
		*pA++ = value;
}


template <class T> void Matrix<T>::setOne()
{
	if ( r!=c )
	{	Fehler(1);
		return;
	}
	
T	*pA = &(*this)[0];
	memset(pA, 0, r*c*sizeof(T));
	for ( INT i=0 ; i<r ; i++ , pA+=r+1)
		*pA = 1;
}



template <class T> 
void Matrix<T>::setOrtho(double phi, INT i1, INT i2)
{
	if ( r!=c )
	{	Fehler(1);
		return;
	}
	
T	*pA = &(*this)[0];
	memset(pA, 0, r*c*sizeof(T));

INT	ii = i1-i2;
	for ( INT i=0 ; i<r ; i++ , pA+=r+1)
		if ( i!=i1 && i!=i2 )
			*pA = 1;
		else
		{
			*pA = (T) cos(phi);
			if ( i==i1 )
				*(pA-ii) =(T) -sin(phi);
			else
				*(pA+ii) =(T) sin(phi);
		}
	
}




template <class T> 
void Matrix<T>::transpose()
{
	if ( r!=c )
	{	Fehler(1);
		return;
	}
	
	for ( INT i=0 ; i<r ; ++i )
		for ( INT j=i+1 ; j<r ; ++j )
		{
		T	h = (*this)(i, j);
			(*this)(i, j) = (*this)(j, i);
			(*this)(j, i) = h;
		}
}



template <> 
void Matrix<double>::inverse()
{
	if ( r!=c )
	{	Fehler(1);
		return;
	}
INT*	ipiv = new INT[r];
INT	info;

	// calculate LU factorization
	FORTRAN_LINKAGE(dgetrf)(&r, &c, p, &r, ipiv, &info);
	assert(!info);
//	cout << "*********** LU info = " << info << " ************" << endl;

INT	lwork = r*r;
double*	work = new double[lwork];
	// calculate inverse
	FORTRAN_LINKAGE(dgetri)(&r, p, &r, ipiv, work, &lwork, &info);
	assert(!info);
//	cout << "*********** INV info = " << info << " ************" << endl;
	transpose();
	delete work;
	delete ipiv;
}





template <class T> Matrix<T> Matrix<T>::operator+()
{
	return *this;
}

template <class T> Matrix<T> Matrix<T>::operator-()
{
INT	i;
Matrix<T>	v(this->getRows(), this->getCols());
	for ( i=0 ; i< (this->getRows())*(this->getCols()) ; i++ )
		v[i] = - (*this)[i];

	return v;
}



template <class T>	void Matrix<T>::Fehler(const INT FehlerNr) const 
{
	puts("Fehler: ");
	switch ( FehlerNr )
	{   case 0	:	puts("not enough memory");
					break;
		case 1	:	puts("Falsche Dimension");
					break;
	}
}


template <class T> Matrix<T> Matrix<T>::Trn()
{
INT	i, j, c;
T	*B;

	c = this->getCols();
	B = new T[(this->getRows())*c];
	for ( i=0 ; i<this->getRows() ; i++ )
	{   *(B+i*(1+c)) = (*this)[i*(1+c)];
		for ( j=i+1 ; j<c ; j++ )
		{	*(B+i*c+j) = (*this)[j*c+i];
			*(B+j*c+i) = (*this)[i*c+j];
		}
	}
	return Matrix<T>(this->getRows(), c, B);
}

template <class T> T Matrix<T>::getTrace()
{
	if ( c!=r )
	{	Fehler(1);
		return 0;
	}

T	t = 0;
	for ( INT i=0 ; i<r ; i++ )
		t += p[i*(1+c)];
	return t;
}


template <class T> double Matrix<T>::getNorm(double Order) const
{
double	N = 0;
T	*pp = p;

	for ( INT i=0 ; i<r*c ; i++ )
		N += pow((double) abs(*pp++), Order);
			
	return pow(N, 1/Order);
}


template <class T> double Matrix<T>::getNorm1() const
{
double	N = 0;
T	*pp = p;

	for ( INT i=0 ; i<r*c ; i++ )
		N += abs(*pp++);
			
	return N;
}


template <class T> double Matrix<T>::getNorm2() const
{
double	N = 0;
T	*pp = p;

	for ( INT i=0 ; i<r*c ; i++ , pp++ )
		N += abs(*pp * *pp);
			
	return sqrt(N);
}


template <class T> double Matrix<T>::getInfNorm() const
{
double	N = 0;
T	*pp = p;
double	n;

	for ( INT i=0 ; i<r*c ; i++ )
		if ( (n=abs(*pp++))>N )
			N = n;
			
	return N;
}




template <class T> Vector<T>
Matrix<T>::GetColVek(const Matrix<T> &A, INT n)
{
INT		i;
	n--;
Vector<T>	v(A.getRows());
	for ( i=0 ; i<A.getRows() ; i++ )
		v[i] = A[n + i*A.getCols()];
    return v;
}

template <class T> Vector<T>
Matrix<T>::GetRowVek(const Matrix<T> &A, INT n)
{
INT		i;
	n--;
Vector<T>	v(A.getRows());
	for ( i=0 ; i<A.getRows() ; i++ )
		v[i] = A[i + n*A.getCols()];
    return v;
}



template <class T>
Matrix<T> Matrix<T>::invert() const
{
Matrix<T>	A(*this);
Matrix<T>	I(r, r);
	I.setOne();

	for ( INT i=1 ; i<r ; i++ )
	{
		if ( A(i, i)==(T)0 )
			continue;
		for ( INT j=i+1 ; j<=r ; j++ )
		{
			if ( A(j, i)==(T)0 )
				continue;
		T	m = -A(j, i)/A(i, i);
			for ( INT k=i ; k<=r ; k++ )
				A(j, k) += A(i, k)*m;
			for ( INT k=1 ; k<=r ; k++ )
				I(j, k) += I(i, k)*m;
		}
	}
	

	for ( INT i=r ; i>1 ; i-- )
	{
		if ( A(i, i)==(T)0 )
			continue;
		for ( INT j=i-1 ; j>=1 ; j-- )
		{
			if ( A(j, i)==(T)0 )
				continue;
		T	m = -A(j, i)/A(i, i);
			A(j, i) += A(i, i)*m;
			for ( INT k=1 ; k<=r ; k++ )
				I(j, k) += I(i, k)*m;
		}
	}

	for ( INT i=1 ; i<=r ; i++ )
	{
	T	m = (T)1/A(i, i);
			for ( INT j=1 ; j<=r ; j++ )
				I(i, j) *= m;
	}

	return I;	
}




template <class T> ostream& operator<<(ostream& s, const Matrix<T> & v)
{
INT	i, j, n;
ios::fmtflags fmt = s.setf(ios::fixed, ios::floatfield);
	s << '[';
	for ( j=n=0 ; j<v.getRows() ; j++ )
	{	s << "[ ";
		for ( i=0 ; i<v.getCols() ; i++ )
			s << setw(12) << setprecision(5) << v[n++] << ' ';
		s << ']';
		if ( j<v.getRows()-1 )
			s << "\n ";
	}
	s << ']';
	s.setf(fmt);
	return s;
}


//*************************************************************************
//*************************************************************************
//*************************************************************************


template class Matrix<INT>;
template class Matrix<float>;
template class Matrix<double>;
template class Matrix<Complex>;

#define INSTANTIATE_TEMPLATE_FUNCTIONS \
template class Matrix<T> operator*(const T, const Matrix<T> &);\
template class Matrix<T> operator+(const Matrix<T> &, const Matrix<T> &);\
template class Matrix<T> operator-(const Matrix<T> &, const Matrix<T> &);\
template class Matrix<T> operator*(const Matrix<T> &, const Matrix<T> &);\
template class Vector<T> operator*(const Matrix<T> &, const Vector<T> &);\
INT operator==(const Matrix<T> &, const Matrix<T> &);\
INT operator!=(const Matrix<T> &, const Matrix<T> &);\
template ostream& operator<<(ostream& s, const Matrix<T> & v);




#define T INT
INSTANTIATE_TEMPLATE_FUNCTIONS
#undef T


#define T float
INSTANTIATE_TEMPLATE_FUNCTIONS
#undef T

#define T double
INSTANTIATE_TEMPLATE_FUNCTIONS
#undef T

#define T Complex
INSTANTIATE_TEMPLATE_FUNCTIONS
#undef T

