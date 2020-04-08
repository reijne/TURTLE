//***********************************************************************
//
//	Name:			BufferedVector.h
//
//	Description:	
//
//	Author:			Michael Hanrath
//
//	Version:		0.0
//
//	Date:			09. Nov 1998
//
//***********************************************************************




/************************************************************************


		Header-File implementiert Vectorrechnung


		V1.0	06.02.93

		V2.0	19.02.96


		(c) by Michael Hanrath


************************************************************************/

#ifndef __BUFFEREDVECTOR_H
#define __BUFFEREDVECTOR_H

#include "../../../config.h"

#include "../Group/Complex.h"
#include <stdarg.h>
#include <stdio.h>
#include <iostream>


/***************************************************************************
****************************************************************************
****************************************************************************/


#include "Vector.h"

class DiskBuffer;


template <class T>
class BufferedVector : private Vector<T> {
public:
	BufferedVector(INT n = 0, INT bufSize = 2*(1 << 20));
	BufferedVector(const DiskBuffer &, INT i, INT bufSize = 2*(1 << 20));
	~BufferedVector();

	BufferedVector(const BufferedVector<T> &);
	BufferedVector & operator = (const BufferedVector<T> &);


	INT	getBufSize() const;

	void	put(const T *);
	void	get(T *) const;


/*	operator BufferedVector<double>();
	operator BufferedVector<Complex>();
*/

	BufferedVector<T> operator * (T) const;

	INT operator==(const BufferedVector<T> &) const;
	INT operator!=(const BufferedVector<T> &) const;


	BufferedVector<T> & operator+=(const BufferedVector<T> &);
	BufferedVector<T> & operator-=(const BufferedVector<T> &);

	BufferedVector<T> & operator*=(T);
	BufferedVector<T> & operator/=(T);

	BufferedVector<T> operator+();
	BufferedVector<T> operator-();



	void	normalize(double Order);
	void	normalize1();
	void	normalize2();
	void	normalizeInf();

	double	getNorm(double Order) const;
	double	getNorm1() const;
	double	getNorm2() const;
	double	getNorm2Sqr() const;
	double	getInfNorm() const;


	T	Sum();

	INT	getDim() const
	{	return	this->N;	}
	

	// for efficiency reasons: non virtual
	const T& operator[](INT n) const;
	// for efficiency reasons: non virtual
	T& operator[](INT n);

	void	clear();


	void	Fehler(INT FehlerNr) const;

protected:
static INT	id;
INT	bufSize;
unsigned INT	shift;
unsigned INT	mask;
char	FileName[1000];
INT	fd;

INT	blockInBuffer;
INT	dirty;

private:
	void	createFile();
	INT	isInBuf(INT i) const;
	void	fillBuf(INT i);
	void	flush();
};


//	friend BufferedVector<Complex> Complex(const BufferedVector<double> &);
//	friend BufferedVector<double> real(const BufferedVector<Complex> &);
//	friend BufferedVector<double> imag(const BufferedVector<Complex> &);

// Wegen besserer Effizienz:
BufferedVector<Complex> operator+(const BufferedVector<Complex> &, const BufferedVector<double> &);
BufferedVector<Complex> operator-(const BufferedVector<Complex> &, const BufferedVector<double> &);
Complex operator*(const BufferedVector<Complex> &, const BufferedVector<double> &);
Complex operator^(const BufferedVector<Complex> &, const BufferedVector<double> &); 		 // diskrete Faltung

BufferedVector<Complex> operator+(const BufferedVector<double> &, const BufferedVector<Complex> &);
BufferedVector<Complex> operator-(const BufferedVector<double> &, const BufferedVector<Complex> &);
Complex operator*(const BufferedVector<double> &, const BufferedVector<Complex> &);
Complex operator^(const BufferedVector<double> &, const BufferedVector<Complex> &); 		 // diskrete Faltung

template <class T>	BufferedVector<T> operator+(const BufferedVector<T> &, const BufferedVector<T> &);
template <class T>	BufferedVector<T> operator-(const BufferedVector<T> &, const BufferedVector<T> &);
template <class T>	T operator*(const BufferedVector<T> &, const BufferedVector<T> &);
double operator*(const BufferedVector<double> &, const BufferedVector<double> &);
float operator*(const BufferedVector<float> &, const BufferedVector<float> &);
template <class T>	T operator^(const BufferedVector<T> &, const BufferedVector<T> &);  		 // diskrete Faltung

template <class T>	BufferedVector<T> operator * (T, const BufferedVector<T> &);

template <class T>	ostream& operator << (ostream& s, const BufferedVector<T> & v);


template <class T>
inline
INT	BufferedVector<T>::getBufSize() const
{	return bufSize;	}


template <class T>
inline
const T& BufferedVector<T>::operator[](INT n) const
{
//	cout << n << " " << (n >> shift) << " " << blockInBuffer << endl;
	if ( !isInBuf(n) )
		((BufferedVector<T> *) this)->fillBuf(n);
	return 	this->p[n & mask];
}

template <class T>
inline
T& BufferedVector<T>::operator[](INT n)
{
	if ( !isInBuf(n) )
		fillBuf(n);
	dirty = 1;
	return 	this->p[n & mask];
}

template <class T>
inline
INT	BufferedVector<T>::isInBuf(INT i) const
{	return (i >> shift) == blockInBuffer;	}


#endif
