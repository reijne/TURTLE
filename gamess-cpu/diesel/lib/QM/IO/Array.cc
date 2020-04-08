//***********************************************************************
//
//	Name:			Array.cc
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


#include "Array.h"

#include <ctype.h>
#include <stdlib.h>

template <class T>
Array<T>::Array(T *_p, T _max)
{
	p = _p;
	max = _max;
}


template <class T>
ostream& operator<<(ostream & s, const Array<T> & a)
{
	for ( INT i=0 ; i<a.getNumber() ; i++ )
		s << a[i] << " ";
	return s;
}

template <class T>
istream& operator>>(istream & s, Array<T> & a)
{
char	buf[1000];
INT	i = 0;
INT	flag;
	do
	{
		while ( s && s.get(buf[i]) && buf[i]!='\n' && i<1000-1 )
			i++;
		if ( !s )
		{
			a.setNumber(0);
			return s;
		}
		buf[i] = 0;
		flag = 0;
		for ( unsigned INT j=0 ; j<strlen(buf) ; j++ )
			if ( buf[j]!=' ' || buf[j]!='\n' )
			{
				flag = 1;
				break;
			}
	} while ( !flag );

INT	n;
char	*p;

	n = i = 0;
	while ( buf[i] )
	{	while ( buf[i] && !isdigit(buf[i]) )
			i++;
		if ( !isdigit(buf[i]) )
			break;
		p = buf+i;
		while ( buf[i] && isdigit(buf[i]) )
			i++;
		if ( buf[i] )
			buf[i++] = 0;
		a[n++] = atoi(p);
		if ( n>=a.getMax() )
			break;		
	}
	a.setNumber(n);
	
	return s;
}



template class Array<INT>;
//template class Array<SHORT_INT>;

template ostream& operator<<(ostream & s, const Array<INT> &);
template istream& operator>>(istream & s, Array<INT> &);

//template ostream& operator<<(ostream & s, const Array<SHORT_INT> &);
//template istream& operator>>(istream & s, Array<SHORT_INT> &);

