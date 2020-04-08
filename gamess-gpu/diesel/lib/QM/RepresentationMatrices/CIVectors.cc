//***********************************************************************
//
//	Name:			CIVectors.cc
//
//	Description:	
//
//					
//
//	Author:			Michael Hanrath
//
//	Version:		0.0
//
//	Date:			10.09.1996
//
//
//
//
//
//***********************************************************************

#include "CIVectors.h"

#include <string>


template <class T>
CIVectors<T>::CIVectors(INT _dim, INT _n, T *_p)
{
	dim = _dim;
	n = _n;
	p = _p;
	shared = NULL;
}


template <class T>
CIVectors<T>::CIVectors(INT _dim, INT _n, SharedMemory::Mode mode)
{
	dim = _dim;
	n = _n;
	shared = new SharedMemory(n*dim*sizeof(T), mode);
	p = (T*) shared->allocate(n*dim, sizeof(T));
}


template <class T>
CIVectors<T>::~CIVectors()
{
	if ( shared )
		delete shared;
}

template <class T>
void	CIVectors<T>::clear()
{
	memset(p, 0, n*dim*sizeof(T));
}

template <class T>
void	CIVectors<T>::setSlave()
{
	if ( shared )
		shared->setSlave();
}

template <class T>
ostream & operator<<(ostream & s, const CIVectors<T> & v)
{
T	*p = v.p;
	for ( INT i=0 ; i<10*0+1*v.dim ; i++ )
	{	for ( INT j=0 ; j<v.n ; j++ )
		{	s << p[j*v.dim];
			if ( j<v.n-1 )
				s << "\t";
		}
		p++;
		s << std::endl;
	}
	return s;
}

template class CIVectors<float>;
template class CIVectors<double>;

template ostream& operator << (ostream& s, const
    CIVectors<double> & v);

template ostream& operator << (ostream& s, const
    CIVectors<float> & v);

