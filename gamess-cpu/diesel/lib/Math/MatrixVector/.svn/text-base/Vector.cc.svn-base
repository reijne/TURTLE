//***********************************************************************
//
//	Name:			Vector.cc
//
//	Description:	
//
//	Author:			Michael Hanrath
//
//	Version:		0.0
//
//	Date:			03. Nov 1998
//
//***********************************************************************


#include "Vector.h"

#include <string>	
#include <stdlib.h>

using namespace std;

template <class T>
Vector<T>::Vector(INT n, T *_p)
{
	N = n;
	p = _p;
}

template <>
Vector<INT>::Vector(INT n, INT a1 ...)
{
INT	i;
va_list	ap;

	N = n;
	p = new INT[n];
	p[0] = a1;

	va_start(ap, a1);
	for ( i=1 ; i<n ; i++ )
		p[i] = va_arg(ap, INT);
	va_end(ap);
}


template <>
Vector<float>::Vector(INT n, float a1 ...)
{
INT	i;
va_list	ap;

	N = n;
	p = new float[n];
	p[0] = a1;

	va_start(ap, a1);
	for ( i=1 ; i<n ; i++ )
		p[i] = va_arg(ap, double);
	va_end(ap);
}

template <>
Vector<double>::Vector(INT n, double a1 ...)
{
INT	i;
va_list	ap;

	N = n;
	p = new double[n];
	p[0] = a1;

	va_start(ap, a1);
	for ( i=1 ; i<n ; i++ )
		p[i] = va_arg(ap, double);
	va_end(ap);
}





/*Vector<INT>::operator Vector<double>()
{
INT	i;
Vector<double>	*A;

	A = new Vector<double>(N);
	for ( i=0 ; i<N ; i++ )
		(*A)[i] = p[i];

	return *A;
}


Vector<float>::operator Vector<double>()
{
INT	i;
Vector<double>	*A;

	A = new Vector<double>(N);
	for ( i=0 ; i<N ; i++ )
		(*A)[i] = p[i];

	return *A;
}


Vector<double>::operator Vector<Complex>()
{
INT	i;
Vector<Complex>	*A;

	A = new Vector<Complex>(N);
	for ( i=0 ; i<N ; i++ )
		(*A)[i] = p[i];

	return *A;
}
*/

/*
Vector<Complex> Complex(Vector<double> & vin)
{
Vector<Complex>	*v;
INT	i;

	v = new Vector<Complex>(vin.getDim());

	for ( i=0 ; i<vin.getDim() ; i++ )
		(*v)[i] = ::Complex(vin[i]);
	return *v;
}
*/

/*
Vector<double> real(const Vector<Complex> & vin)
{
Vector<double>	*v;
INT	i;

	v = new Vector<double>(vin.getDim());

	for ( i=0 ; i<vin.getDim() ; i++ )
		(*v)[i] = ::real(vin[i]);
	return *v;
}

Vector<double> imag(const Vector<Complex> & vin)
{
Vector<double>	*v;
INT	i;

	v = new Vector<double>(vin.getDim());

	for ( i=0 ; i<vin.getDim() ; i++ )
		(*v)[i] = ::imag(vin[i]);
	return *v;
}


*/

//---------------------------------------------------------------

Vector<Complex> operator+(const Vector<Complex> &v1, const Vector<double> &v2)
{
INT	i;
	if ( v1.getDim() != v2.getDim() )
	{	v1.Fehler(1);
Vector<Complex> 	v;
		return v;
	}
Vector<Complex>	v(v1.getDim());
	for ( i=0 ; i<v1.getDim() ; i++ )
		v[i] = v1[i] + v2[i];
	return v;
}

Vector<Complex> operator-(const Vector<Complex> &v1, const Vector<double> &v2)
{
INT	i;
	if ( v1.getDim() != v2.getDim() )
	{	v1.Fehler(1);
Vector<Complex> 	v;
		return v;
	}
Vector<Complex>	v(v1.getDim());
	for ( i=0 ; i<v1.getDim() ; i++ )
		v[i] = v1[i] - v2[i];
	return v;
}

Complex operator*(const Vector<Complex> &v1, const Vector<double> &v2)
{
INT	i;
Complex	sum = 0;
	if ( v1.getDim() != v2.getDim() )
	{	v1.Fehler(1);
		return sum;
	}
	for ( i=0 ; i<v1.getDim() ; i++ )
		sum += v1[i] * v2[i];
	return sum;
}

Complex operator^(const Vector<Complex> &v1, const Vector<double> &v2)
{
INT	i, n;
Complex	sum = 0;
	if ( v1.getDim() != v2.getDim() )
	{	v1.Fehler(1);
		return sum;
	}
	n = v1.getDim();
	for ( i=0 ; i<n ; i++ )
		sum += v1[n-i-1] * v2[i];
	return sum;
}

Vector<Complex> operator%(const Vector<Complex> &v1, const Vector<double> &v2)
{
	if ( v1.getDim() != v2.getDim() )
	{	v1.Fehler(1);
Vector<Complex> 	v;
		return v;
	}
	if ( v1.getDim() != 3 )
	{	v1.Fehler(2);
Vector<Complex> 	v;
		return v;
	}
Vector<Complex>	v(3);
	v[0] =   v1[1] * v2[2]  -  v1[2] * v2[1];
	v[1] = - v1[0] * v2[2]  +  v1[2] * v2[0];
	v[2] =   v1[0] * v2[1]  -  v1[1] * v2[0];
	return v;
}



Vector<Complex> operator+(const Vector<double> &v1, const Vector<Complex> &v2)
{	return	v2+v1;	};

Vector<Complex> operator-(const Vector<double> &v1, const Vector<Complex> &v2)
{	return	-(v2-v1);	};

Complex operator*(const Vector<double> &v1, const Vector<Complex> &v2)
{	return	v2*v1;	};

Complex operator^(const Vector<double> &v1, const Vector<Complex> &v2)			// diskrete Faltung
{	return	v2^v1;	};

Vector<Complex> operator%(const Vector<double> &v1, const Vector<Complex> &v2)	// Vectormultiplikation
{	return -(v2%v1);	};



//---------------------------------------------------------------






//*********************************************************************
//*********************************************************************
//*********************************************************************
//*********************************************************************
//*********************************************************************
//*********************************************************************

template <class T> Vector<T>::Vector(INT n)
{
	N = n;
	p = new T[n];
}




template <class T> Vector<T>::~Vector()
{
	N = 0;
	if ( p )
		delete[] p;
}




template <class T> Vector<T>& Vector<T>::operator=(const Vector<T> & v1)
{
INT	i;

	if ( v1.getDim() != N )
	{	if ( N )
			delete	p;
		p = new T[v1.getDim()];
		N = v1.getDim();
	}
	for ( i=0 ; i<v1.getDim() ; i++ )
		p[i] = v1[i];
	return *this;
}


template <class T> Vector<T>::Vector(const Vector<T> & v)
{
	N = v.getDim();
	p = new T[N];
	memcpy(p, v.getP(), N*sizeof(T));
}



template <>
Vector<float>::Vector(INT n, INT a1 ...)
{
INT	i;
va_list	ap;

	N = n;
	p = new float[n];
	p[0] = a1;

	va_start(ap, a1);
	for ( i=1 ; i<n ; i++ )
		p[i] = va_arg(ap, INT);
	va_end(ap);
}


template <>
Vector<double>::Vector(INT n, INT a1 ...)
{
INT	i;
va_list	ap;

	N = n;
	p = new double[n];
	p[0] = a1;

	va_start(ap, a1);
	for ( i=1 ; i<n ; i++ )
		p[i] = va_arg(ap, INT);
	va_end(ap);
}


template <class Complex>
Vector<Complex>::Vector(INT n, INT a1 ...)
{
INT	i;
va_list	ap;

	N = n;
	p = new Complex[n];
	p[0] = a1;

	va_start(ap, a1);
	for ( i=1 ; i<n ; i++ )
		p[i] = va_arg(ap, INT);
	va_end(ap);
}



template <class T> Vector<T> operator*(T a, const Vector<T> &v1)
{
INT	i;
Vector<T>	v(v1.getDim());
	for ( i=0 ; i<v1.getDim() ; i++ )
		v[i] = a*v1[i];
	return v;
}


template <class T> 
Vector<T> Vector<T>::operator*(T a) const
{
INT	i;
Vector<T>	v(getDim());
	for ( i=0 ; i<N ; i++ )
		v[i] = a*p[i];
	return v;
}


template <class T> Vector<T> operator+(const Vector<T> &v1, const Vector<T> &v2)
{
INT	i;
	if ( v1.getDim() != v2.getDim() )
	{	v1.Fehler(1);
Vector<T> 	v;
		return v;
	}
Vector<T>	v(v1.getDim());
	for ( i=0 ; i<v1.getDim() ; i++ )
		v[i] = v1[i] + v2[i];
	return v;
}

template <class T> Vector<T> operator-(const Vector<T> &v1, const Vector<T> &v2)
{
INT	i;
	if ( v1.getDim() != v2.getDim() )
	{	v1.Fehler(1);
Vector<T> 	v;
		return v;
	}
Vector<T>	v(v1.getDim());
	for ( i=0 ; i<v1.getDim() ; i++ )
		v[i] = v1[i] - v2[i];
	return v;
}

template <class T> T operator*(const Vector<T> &v1, const Vector<T> &v2)
{
INT     i;
T       sum = 0;
        if ( v1.getDim() != v2.getDim() )
        {       v1.Fehler(1);
                return sum;
        }
        for ( i=0 ; i<v1.getDim() ; i++ )
                sum += v1[i] * v2[i];
        return sum;
}

double operator*(const Vector<double> &v1, const Vector<double> &v2)
{
INT	i;
double	sum = 0;
	if ( v1.getDim() != v2.getDim() )
	{	v1.Fehler(1);
		return sum;
	}
	for ( i=0 ; i<v1.getDim() ; i++ )
		sum += v1[i] * v2[i];
	return sum;
}

float operator*(const Vector<float> &v1, const Vector<float> &v2)
{
INT	i;
double	sum = 0;
	if ( v1.getDim() != v2.getDim() )
	{	v1.Fehler(1);
		return sum;
	}
	for ( i=0 ; i<v1.getDim() ; i++ )
		sum += v1[i] * v2[i];
	return sum;
}

template <class T> T operator^(const Vector<T> &v1, const Vector<T> &v2)
{
INT	i, n;
T 	sum = 0;
	if ( v1.getDim() != v2.getDim() )
	{	v1.Fehler(1);
		return sum;
	}
	n = v1.getDim();
	for ( i=0 ; i<n ; i++ )
		sum += v1[n-i-1] * v2[i];
	return sum;
}

template <class T> Vector<T> operator%(const Vector<T> &v1, const Vector<T> &v2)
{
	if ( v1.getDim() != v2.getDim() )
	{	v1.Fehler(1);
Vector<T> 	v;
		return v;
	}
	if ( v1.getDim() != 3 )
	{	v1.Fehler(2);
Vector<T> 	v;
		return v;
	}
Vector<T>	v(3);
	v[0] =   v1[1] * v2[2]  -  v1[2] * v2[1];
	v[1] = - v1[0] * v2[2]  +  v1[2] * v2[0];
	v[2] =   v1[0] * v2[1]  -  v1[1] * v2[0];
	return v;
}



template <class T>	T Vector<T>::Sum()
{
INT	i;
T	sum = 0;
	for ( i=0 ; i<N ; i++ )
		sum += (*this)[i];
	return sum;
}



template <class T> 	void Vector<T>::Fehler(INT FehlerNr) const
{
	puts("Fehler: ");
	switch ( FehlerNr )
	{	case 1	:	puts("Falsche Dimension");
					break;
		case 2	:	puts("Falsche Dimension (!=3) bei Vectorprodukt");
					break;
	}
}


template <class T> 
INT Vector<T>::operator==(const Vector<T> & v) const
{
INT	i;
	if ( getDim() != v.getDim() )
	{	Fehler(1);
		return 0;
	}
	for ( i=0 ; i<getDim() ; i++ )
		if ( (*this)[i] != v[i] )
			return 0;
	return 1;
}


template <class T> 
INT Vector<T>::operator!=(const Vector<T> & v) const
{
INT	i;
	if ( getDim() != v.getDim() )
	{	Fehler(1);
		return 0;
	}
	for ( i=0 ; i<getDim() ; i++ )
		if ( (*this)[i] != v[i] )
			return 1;
	return 0;
}



template <class T> Vector<T> & Vector<T>::operator+=(const Vector<T> &v1)
{
INT	i;
	if ( this->getDim() != v1.getDim() )
	{	v1.Fehler(1);
		return *this;
	}
	for ( i=0 ; i<v1.getDim() ; i++ )
		(*this)[i] += v1[i];
	return *this;
}

template <class T> Vector<T> & Vector<T>::operator-=(const Vector<T> &v1)
{
INT	i;
	if ( this->getDim() != v1.getDim() )
	{	v1.Fehler(1);
		return *this;
	}
	for ( i=0 ; i<v1.getDim() ; i++ )
		(*this)[i] -= v1[i];
	return *this;
}

template <class T> Vector<T> & Vector<T>::operator%=(const Vector<T> &v2)
{
Vector<T>	v1(3);
	if ( v1.getDim() != v2.getDim() )
	{	v1.Fehler(1);
		return *this;
	}
	if ( v1.getDim() != 3 )
	{	v1.Fehler(2);
		return *this;
	}
	v1 = *this;
	(*this)[0] =   v1[1] * v2[2]  -  v1[2] * v2[1];
	(*this)[1] = - v1[0] * v2[2]  +  v1[2] * v2[0];
	(*this)[2] =   v1[0] * v2[1]  -  v1[1] * v2[0];
	return *this;
}



template <class T> Vector<T> Vector<T>::operator+()
{
	return *this;
}

template <class T> Vector<T> Vector<T>::operator-()
{
INT	i;
Vector<T>	v(this->getDim());
	for ( i=0 ; i< this->getDim() ; i++ )
		v[i] = - (*this)[i];

	return v;
}



template <class T> Vector<T> & Vector<T>::operator*=(T a)
{
	for ( INT i=0 ; i<getDim() ; i++ )
		p[i] *= a;
	return *this;
}

template <class T> Vector<T> & Vector<T>::operator/=(T a)
{
	for ( INT i=0 ; i<getDim() ; i++ )
		p[i] /= a;
	return *this;
}


template <class T> void Vector<T>::normalize(double Order)
{
	*this /= (T)getNorm(Order);
}

template <class T> void Vector<T>::normalize1()
{
	*this /= (T)getNorm1();
}

template <class T> void Vector<T>::normalize2()
{
	*this /= (T)getNorm2();
}

template <class T> void Vector<T>::normalizeInf()
{
	*this /= (T)getInfNorm();
}



template <class T> double Vector<T>::getNorm(double Order) const
{
double	Norm = 0;
T	*pp = p;

	for ( INT i=0 ; i<N ; i++ )
		Norm += pow((double)abs(*pp++), Order);
			
	return pow(Norm, 1/Order);
}


	


template <class T> double Vector<T>::getNorm1() const
{
double	Norm = 0;
T	*pp = p;

	for ( INT i=0 ; i<N ; i++ )
		Norm += abs(*pp++);
			
	return Norm;
}


template <class T> double Vector<T>::getNorm2() const
{
double	Norm = 0;
T	*pp = p;

	for ( INT i=0 ; i<N ; i++ , pp++ )
		Norm += abs(*pp * *pp);
			
	return sqrt(Norm);
}


template <class T> double Vector<T>::getNorm2Sqr() const
{
double	Norm = 0;
T	*pp = p;

	for ( INT i=0 ; i<N ; i++ , pp++ )
		Norm += abs(*pp * *pp);
			
	return Norm;
}


template <class T> double Vector<T>::getInfNorm() const
{
double	Norm = 0;
T	*pp = p;
double	n;

	for ( INT i=0 ; i<N ; i++ )
		if ( (n=abs(*pp++))>Norm )
			Norm = n;
			
	return Norm;
}


template <class T>
void	Vector<T>::clear()
{
	memset(p, 0, N*sizeof(T));
}


template <class T> ostream& operator<<(ostream& s, const Vector<T> & v)
{
INT	i;
	s << "[ ";
	for ( i=0 ; i<v.getDim() ; i++ )
	{	s << v[i] << ' ';
	}
	s << ']';
	return s;
}

template <class T> void Vector<T>::setDim(INT n)
{
        N = n;
        p = new T[n];
}


//*************************************************************************
//*************************************************************************
//*************************************************************************
//*************************************************************************


template class Vector<INT>;
template class Vector<float>;
template class Vector<double>;
template class Vector<Complex>;


#define INSTANTIATE_TEMPLATE_FUNCTIONS \
template class Vector<T> operator*(T, const Vector<T> &);\
template class Vector<T> operator+(const Vector<T> &, const Vector<T> &);\
template class Vector<T> operator-(const Vector<T> &, const Vector<T> &);\
template class Vector<T> operator%(const Vector<T> &, const Vector<T> &);\
template ostream& operator << (ostream& s, const Vector<T> & v);\
T operator*(const Vector<T> &, const Vector<T> &);\
T operator^(const Vector<T> &, const Vector<T> &);


/*
// Macke im GNU C++ 2.7.0:
// Instantiierung der Nicht-Memberfunctions, die integralen Typ
// zurückgeben, funktioniert nicht! Daher folgende
// "Zwangs-Instantiierung":

static void	dummy()
{
Complex	h;
Vector<INT>	vi;
	vi!=vi;
	vi==vi;
	h = vi*vi;
	h = vi^vi;
	
Vector<float>	vf;
	vf!=vf;
	vf==vf;
	h = vf*vf;
	h = vf^vf;
	
Vector<double>	vd;
	vd!=vd;
	vd==vd;
	h = vd*vd;
	h = vd^vd;
	
Vector<Complex>	vc;
	vc!=vc;
	vc==vc;
	h = vc*vc;
	h = vc^vc;
	
}
*/


//template class INT operator!=(const Vector<INT> &, const Vector<INT> &);

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

