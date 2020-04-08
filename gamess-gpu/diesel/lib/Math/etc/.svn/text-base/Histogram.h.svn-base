//***********************************************************************
//
//	Name:			Histogram.h
//
//	Description:	
//
//	Author:			Michael Hanrath
//
//	Version:		0.0
//
//	Date:			22.04.1997
//
//
//
//
//***********************************************************************

#ifndef __Histogram_h
#define __Histogram_h

#include "../../../config.h"

#include <iostream>
using std::ostream;

#include <math.h>

template <class T>
class Histogram {
public:
enum ScaleMode { Linear, Logarithmic };
	Histogram(
		T minValue,
		T maxValue, 
//	GNU C++ 2.7.2 work around
//		ScaleMode mode = Linear,
		INT mode = 0,
		INT steps = 40);
	~Histogram();
	
	
	Histogram(const Histogram<T> &);
	Histogram<T> & operator = (const Histogram<T> &);
	
	
	//------------------------------------------------------------------

	INT	count(T value);
	
	//------------------------------------------------------------------
	
	INT	getSteps() const;
	
	INT	getFrequencyByValue(T value) const;
	INT getFrequencyByIntervalNo(INT i) const;
	
	INT	getIndex(T value) const;

	INT getNumberOfCounts() const;
	INT	getNumberOfHitsBelow() const;
	INT	getNumberOfHitsAbove() const;
	
	T getIntervalMin(INT i) const;
	T getIntervalMax(INT i) const;
	
	T getMinValue() const;
	T getMaxValue() const;
	
	//------------------------------------------------------------------

	void	clear();
	
	//------------------------------------------------------------------
	
	void	printPrettyASCII(INT cols = 80);
        template <class TT>
	friend ostream& operator<< (ostream& s, const Histogram<TT>);
	
	//------------------------------------------------------------------
	
	
private:
	
INT	steps;						// number of intervals
INT	counts;						// counts
T	minValue;					// minimum value mapped
T	maxValue;					// maximum value mapped
unsigned INT	*data;			// interval counters
INT	below;						// hits below minimum value
INT	above;						// hits above maximum value
double	scale;					// scaling factor for index calculation
ScaleMode	mode;				// mode of scaling
};



template <class T>
inline
INT	Histogram<T>::getIndex(T value) const
{
	switch ( mode ) {
	case Linear:
		return (INT) ((value-minValue)*scale + 0.5);

	case Logarithmic:
		return (INT) (value<=0 ? 0 : (log(fabs((double)value)/minValue)*scale));
	}
	return 0;
}


template <class T>
inline
INT	Histogram<T>::count(T value)
{
	counts++;
INT	i = getIndex(value);
	if ( i<0 )
	{
		below++;
		return counts;
	}
	if ( i>=steps )
	{
		above++;
		return counts;
	}
	data[i]++;
	return counts;
}


template <class T>
inline
INT	Histogram<T>::getSteps() const
{	return	steps;	}

template <class T>
inline
INT	Histogram<T>::getNumberOfCounts() const
{	return  counts;	}

template <class T>
inline
INT	Histogram<T>::getNumberOfHitsBelow() const
{	return  below;	}


template <class T>
inline
INT	Histogram<T>::getNumberOfHitsAbove() const
{	return  above;	}


template <class T>
inline
INT	Histogram<T>::getFrequencyByValue(T value) const
{
INT	i = getIndex(value);
	if ( i<0 )
		return below;
	if ( i>steps )
		return above;
	return data[i];
}

template <class T>
inline
INT	Histogram<T>::getFrequencyByIntervalNo(INT i) const
{	return  data[i];	}



template <class T>
inline
T	Histogram<T>::getIntervalMin(INT i) const
{
	switch ( mode ) {
	case Linear:
		return  (T) (i/scale+minValue);

	case Logarithmic:
		return  (T) (exp(i/scale)*minValue);
	}
	return 0;
}
	


template <class T>
inline
T	Histogram<T>::getIntervalMax(INT i) const
{
	switch ( mode ) {
	case Linear:
		return  (T) ((i+1)/scale+minValue);

	case Logarithmic:
		return  (T) (exp((i+1)/scale)*minValue);
	}
	return 0;
}



template <class T>
inline
T	Histogram<T>::getMinValue() const
{	return  minValue;	}


template <class T>
inline
T	Histogram<T>::getMaxValue() const
{	return  maxValue;	}



#endif


