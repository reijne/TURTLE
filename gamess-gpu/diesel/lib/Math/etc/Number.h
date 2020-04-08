//***********************************************************************
//
//	Name:			Number.h
//
//	Description:	
//
//	Author:			Michael Hanrath
//
//	Version:		0.0
//
//	Date:			21. June 1999
//
//***********************************************************************

#ifndef __Number_H
#define __Number_H

#include "../../../config.h"

#include <iostream>

template <class T>
class Number {
public:
	Number() : v(0), init(0) {}
	Number(const Number &a) : v(a.v), init(a.init) {}
	Number(T a, INT i = 1) : v(a), init(i) {}
//	operator T() const	{	return v;	}

	T	value() const	{	return v;	}

	INT	isNumber() const {	return init;	}
	void setNumber(INT i) {	init = i;	}
	
	
	
	Number operator + (Number a) const
	{
	Number	n(*this);
		n.init &= a.init;
		n.v += a.v;
		return n;
	}
			
	Number operator - (Number a) const
	{
	Number	n(*this);
		n.init &= a.init;
		n.v -= a.v;
		return n;
	}
			
	Number operator * (Number a) const
	{
	Number	n(*this);
		n.init &= a.init;
		n.v *= a.v;
		return n;
	}
			
	Number operator / (Number a) const
	{
	Number	n(*this);
		n.init &= a.init;
		if ( !n.init )
			return n;
 		n.v /= a.v;
		return n;
	}
			
	Number operator - () const
	{
		return Number(-v, init);
	}

	Number invert () const
	{
		return N(1/v, init);
	}

	INT operator == (Number a) const
	{		return init && a.init && v==a.v;	}
			
	INT operator != (Number a) const
	{		return !init || !a.init || v!=a.v;	}
			
	
	
	Number & operator += (Number a)
	{
		init &= a.init;
		v += a.v;
		return *this;
	}

	Number & operator -= (Number a)
	{
		init &= a.init;
		v -= a.v;
		return *this;
	}

	Number & operator *= (Number a)
	{
		init &= a.init;
		v *= a.v;
		return *this;
	}

	Number & operator /= (Number a)
	{
		init &= a.init;
		if ( !this->n.init )
			return *this;
		v /= a.v;
		return *this;
	}

	
	template <class TT> friend Number<TT> operator + (TT, const Number<TT> &);
	template <class TT> friend Number<TT> operator - (TT, const Number<TT> &);
	template <class TT> friend Number<TT> operator * (TT, const Number<TT> &);
	template <class TT> friend Number<TT> operator / (TT, const Number<TT> &);

	template <class TT> friend Number<TT> & operator += (TT, const Number<TT> &);
	template <class TT> friend Number<TT> & operator -= (TT, const Number<TT> &);
	template <class TT> friend Number<TT> & operator *= (TT, const Number<TT> &);
	template <class TT> friend Number<TT> & operator /= (TT, const Number<TT> &);


	template <class TT> friend ostream& operator << (ostream& s, const Number<TT> &);

//	Number & operator = (Number a)	{	v = a;	init = 1;	return *this;	}


private:
T	v;
INT	init;
};


template <class T>
inline
ostream& operator<<(ostream& s, const Number<T> & n)
{
	if ( n.init )
		s << n.v;
	else
		s << "---";
	return s;
}

template <class T>
inline
Number<T> operator + (T a, const Number<T> & b)
{
	return 	b + a;
}

template <class T>
inline
Number<T> operator - (T a, const Number<T> & b)
{
	return 	-(b - a);
}

template <class T>
inline
Number<T> operator * (T a, const Number<T> & b)
{
	return 	b * a;
}

template <class T>
inline
Number<T> operator / (T a, const Number<T> & b)
{
	return 	invert(b / a);
}



#endif
