//***********************************************************************
//
//	Name:			Vector.h
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




/************************************************************************


		Header-File implementiert Vectorrechnung


		V1.0	06.02.93

		V2.0	19.02.96


		(c) by Michael Hanrath


************************************************************************/

#ifndef __VECTOR_H
#define __VECTOR_H

#include "../../../config.h"

#include "../Group/Complex.h"
#include <stdarg.h>
#include <stdio.h>
#include <iostream>
using std::ostream;

/***************************************************************************
****************************************************************************
****************************************************************************/

template <class T>
class Vector {
public:
	Vector(INT n = 0);
	Vector(INT n, T*);
	Vector(INT n, INT a1 ...);
	Vector(INT n, float a1 ...);
	Vector(INT n, double a1 ...);

	Vector& operator=(const Vector<T> &);
	Vector(const Vector<T> &);



	~Vector();


/*	operator Vector<double>();
	operator Vector<Complex>();
*/


	Vector<T> operator * (T) const;

	INT operator==(const Vector<T> &) const;
	INT operator!=(const Vector<T> &) const;


	Vector<T> & operator+=(const Vector<T> &);
	Vector<T> & operator-=(const Vector<T> &);
	Vector<T> & operator%=(const Vector<T> &);

	Vector<T> & operator*=(T);
	Vector<T> & operator/=(T);

	Vector<T> operator+();
	Vector<T> operator-();



	void	normalize(double Order);
	void	normalize1();
	void	normalize2();
	void	normalizeInf();

	double	getNorm(double Order) const;
	double	getNorm1() const;
	double	getNorm2() const;
	double	getNorm2Sqr() const;
	double	getInfNorm() const;

	void	setDim(INT n);


	T	Sum();

	INT	getDim() const
	{	return	N;	}

	T * getP() const
	{	return p;	}

	T& operator[](const INT n)
 	{	return 	p[n];	}

	const T& operator[](const INT n) const
 	{	return 	p[n];	}

	void	clear();


	void	Fehler(INT FehlerNr) const;

protected:
	INT		N;
	T		*p;


};


//	friend Vector<Complex> Complex(const Vector<double> &);
//	friend Vector<double> real(const Vector<Complex> &);
//	friend Vector<double> imag(const Vector<Complex> &);

// Wegen besserer Effizienz:
Vector<Complex> operator+(const Vector<Complex> &, const Vector<double> &);
Vector<Complex> operator-(const Vector<Complex> &, const Vector<double> &);
Complex operator*(const Vector<Complex> &, const Vector<double> &);
Complex operator^(const Vector<Complex> &, const Vector<double> &); 		 // diskrete Faltung
Vector<Complex> operator%(const Vector<Complex> &, const Vector<double> &);  // Vectormultiplikation

Vector<Complex> operator+(const Vector<double> &, const Vector<Complex> &);
Vector<Complex> operator-(const Vector<double> &, const Vector<Complex> &);
Complex operator*(const Vector<double> &, const Vector<Complex> &);
Complex operator^(const Vector<double> &, const Vector<Complex> &); 		 // diskrete Faltung
Vector<Complex> operator%(const Vector<double> &, const Vector<Complex> &);  // Vectormultiplikation

template <class T>	Vector<T> operator+(const Vector<T> &, const Vector<T> &);
template <class T>	Vector<T> operator-(const Vector<T> &, const Vector<T> &);
template <class T>	T operator*(const Vector<T> &, const Vector<T> &);
double operator*(const Vector<double> &, const Vector<double> &);
float operator*(const Vector<float> &, const Vector<float> &);
template <class T>	T operator^(const Vector<T> &, const Vector<T> &);  		 // diskrete Faltung
template <class T>	Vector<T> operator%(const Vector<T> &, const Vector<T> &);   // Vectormultiplikation

template <class T>	Vector<T> operator * (T, const Vector<T> &);

template <class T>	ostream& operator << (ostream& s, const Vector<T> & v);

#endif
