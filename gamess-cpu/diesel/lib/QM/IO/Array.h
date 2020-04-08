//***********************************************************************
//
//	Name:			Array.h
//
//	Description:	array stream IO
//
//	Author:			Michael Hanrath
//
//	Version:		0.0
//
//	Date:			12.09.1996
//
//
//
//
//
//***********************************************************************

#ifndef __ARRAY_H
#define __ARRAY_H

#include "../../../config.h"

#include <iostream>
using std::ostream;
using std::istream;

template <class T> class Array;
template <class T> ostream& operator<< (ostream & s, const Array<T> &);
template <class T> istream& operator>> (istream & s, Array<T> &);

template <class T>
class Array {
public:
	Array(T *p, T max);
	~Array() {}

	//----------------------------------------------------------------

	T & operator [] (INT i) const;
	
	//----------------------------------------------------------------

	INT getNumber() const;
	void	setNumber(INT n);
	
	T	getMax() const;

	//----------------------------------------------------------------

	friend	ostream& operator<< <T> (ostream & s, const Array<T> &);
	friend	istream& operator>> <T> (istream & s, Array<T> &);
	
private:
T	max;
INT	n;
T	*p;	
};

template <class T>	ostream& operator<< (ostream & s, const Array<T> &);
template <class T>	istream& operator>> (istream & s, Array<T> &);



template <class T>
inline
T & Array<T>::operator [] (INT i) const
{	return p[i];	}

template <class T>
inline
INT	Array<T>::getNumber() const
{	return n;	}

template <class T>
inline
void	Array<T>::setNumber(INT _n)
{	n = _n;	}

template <class T>
inline
T	Array<T>::getMax() const
{	return max;	}


#endif

