//***********************************************************************
//
//	Name:			Histogram.cc
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

#include "Histogram.h"

#include <string>
#include <iomanip>

using namespace std;

template <class T>
Histogram<T>::Histogram(
	T _minValue, T _maxValue,
//	GNU C++ 2.7.2 work around
//	Histogram<T>::ScaleMode _mode,
	INT _mode,
	INT _steps)
{
	minValue = _minValue;
	maxValue = _maxValue;
	steps = _steps;
//	GNU C++ 2.7.2 work around
	mode = (ScaleMode) _mode;
	counts = below = above = 0;
	data = new unsigned INT[steps];
	memset(data, 0, steps*sizeof(unsigned INT));
	switch ( mode ) {
	case Linear:
		scale = steps/(maxValue-minValue);
		break;

	case Logarithmic:
		scale = steps/log((double)maxValue/(double)minValue);
		break;
	}
}


template <class T>
Histogram<T>::~Histogram()
{
	delete data;
}


template <class T>
void	Histogram<T>::clear()
{
	counts = below = above = 0;
	memset(data, 0, steps*sizeof(unsigned INT));
}

template <class T>
void	Histogram<T>::printPrettyASCII(INT cols)
{
const INT	leftSpace = 22;
const INT	rightSpace = 10;

unsigned INT	maxV = 0;
	for ( INT i=0 ; i<steps ; i++ )
		if ( maxV<data[i] )
			maxV = data[i];
			
double	scale = (1.0*(cols - (leftSpace+rightSpace)))/maxV;

INT	prec = 0;
ios::fmtflags format  = ios::fixed;
	switch ( mode ) {
	case Linear:
		format = ios::fixed;
		prec = 4;
		break;

	case Logarithmic:
		format = ios::scientific;
		prec = 3;
		break;

	default:
		format = ios::fixed;
		prec = 4;
		break;
	}

	cout << counts << " counts:" << endl;
	cout << "         >= " << setw(8) << setiosflags(format)
		<< setprecision(prec) << maxValue << " : " 
		<< setw(8) << above << endl;
	for ( INT i=steps-1 ; i>=0 ; i-- )
	{
		cout << "[" << setw(8) << setiosflags(format) << setprecision(prec) 
			<< getIntervalMin(i) << " - " 
			<< setw(8) << setiosflags(format) << setprecision(prec) 
			<< getIntervalMax(i) << "): ";
	INT	j;
		for ( j=0 ; j<(INT) (data[i]*scale) ; j++ )
			cout << "=";
		for ( ; j<cols - (leftSpace+rightSpace) ; j++ )
			cout << " ";
		cout << " " << setw(8) << getFrequencyByIntervalNo(i) << endl;
	}
	cout << "          < " << setw(8)
		<< setprecision(prec)<< setiosflags(format) << minValue << " : " 
		<< setw(8) << below << endl;
}


template <class T>
ostream& operator<<(ostream& s, const Histogram<T>)
{
}



template class Histogram<double>;
template class Histogram<INT>;
