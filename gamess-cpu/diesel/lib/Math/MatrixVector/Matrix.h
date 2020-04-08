//***********************************************************************
//
//	Name:			Matrix.h
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




/************************************************************************


		Header-File implementiert Matrizenrechnung


		V1.0	16.05.93

		V2.0	19.02.96


		(c) by Michael Hanrath


************************************************************************/

#ifndef __MATRIX_H
#define __MATRIX_H

#include "../../../config.h"


#include "../Group/Complex.h"
#include <stdarg.h>
#include <stdio.h>

#include "Vector.h"

/***************************************************************************
****************************************************************************
****************************************************************************/

template <class T>
class Matrix {
public:
	Matrix(INT r = 0, INT c = 0);
	Matrix(INT r, INT c, const T *p);
	Matrix(INT r, INT c, INT a1 ...);
	Matrix(INT r, INT c, float a1 ...);
	Matrix(INT r, INT c, double a1 ...);
	virtual ~Matrix();

/*	operator Matrix<double>();
	operator Matrix<Complex>();
*/

	Matrix<INT> integer(const Matrix<double> &);




//	Vector<T> operator*(Vector<T> &, Matrix<T> &);



	Matrix<T> & operator+=(const Matrix<T> &);
	Matrix<T> & operator-=(const Matrix<T> &);
	Matrix<T> & operator*=(const Matrix<T> &);

	Matrix<T> & operator*=(const T &);


	Matrix<T> invert() const;


	Matrix<T> operator+();
	Matrix<T> operator-();

	Matrix<T> Diagon(const Matrix<T> &);

	Matrix& operator=(const Matrix<T> &);
	Matrix(const Matrix<T> &);




	Matrix Trn();
	T	getTrace();


	Vector<T> GetColVek(const Matrix<T> &, INT n);
	Vector<T> GetRowVek(const Matrix<T> &, INT n);


	double	getNorm(double Order) const;
	double	getNorm1() const;
	double	getNorm2() const;
	double	getInfNorm() const;


	void	setFill(const T & value);
	void	setOne();
	void	setOrtho(double phi, INT i1, INT i2);

	void	transpose();
	void	inverse();

	INT	getCols() const 
	{	return	c;	}

	INT	getRows() const 
	{	return	r;	}

	T * getP() const
	{	return p;	}

	virtual T& operator[](const INT n) const 
 	{	return 	p[n];	}

	T& operator()(const INT i, const INT j) const 
	{	return 	p[i*c + j];	}

	void	Fehler(const INT FehlerNr) const ;

protected:
	INT		r, c;
	T		*p;
};


//	friend Matrix<Complex> Complex(const Matrix<double> &);
//	friend Matrix<double> real(const Matrix<Complex> &);
//	friend Matrix<double> imag(const Matrix<Complex> &);


Matrix<Complex> operator+(const Matrix<Complex> &, const Matrix<double> &);
Matrix<Complex> operator-(const Matrix<Complex> &, const Matrix<double> &);
Matrix<Complex> operator*(const Matrix<Complex> &, const Matrix<double> &);

Matrix<Complex> operator+(const Matrix<double> &, const Matrix<Complex> &);
Matrix<Complex> operator-(const Matrix<double> &, const Matrix<Complex> &);
Matrix<Complex> operator*(const Matrix<double> &, const Matrix<Complex> &);

Vector<Complex> operator*(const Matrix<double> &, const Vector<Complex> &);
Vector<Complex> operator*(const Matrix<Complex> &, const Vector<double> &);

template <class T>	Vector<T> operator*(const Matrix<T> &, const Vector<T> &);
template <class T>	Matrix<T> operator*(const T, const Matrix<T> &);

template <class T>	Matrix<T> operator+(const Matrix<T> &, const Matrix<T> &);
template <class T>	Matrix<T> operator-(const Matrix<T> &, const Matrix<T> &);
template <class T>	Matrix<T> operator*(const Matrix<T> &, const Matrix<T> &);

template <class T>	INT operator==(const Matrix<T> &, const Matrix<T> &);
template <class T>	INT operator!=(const Matrix<T> &, const Matrix<T> &);

template <class T>	ostream& operator<<(ostream& s, const Matrix<T> & v);

#endif
