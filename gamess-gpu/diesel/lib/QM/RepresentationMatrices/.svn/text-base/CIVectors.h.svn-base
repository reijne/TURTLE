//***********************************************************************
//
//	Name:			CIVectors.h
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

#ifndef __CIVECTORS_H
#define __CIVECTORS_H

#include "../../../config.h"

#include <iostream>
using std::ostream;

#include "../Parallel/SharedMemory.h"

template <class T> class CIVectors;
template <class T> ostream & operator<< (ostream & s, const CIVectors<T> &);

template <class T>
class CIVectors {
public:
	CIVectors(INT dim, INT n, T *p);
	CIVectors(INT dim, INT n, SharedMemory::Mode mode = SharedMemory::StandAlone);
	~CIVectors();
	
	
	T &	operator() (INT dim, INT n);
	T &	operator() (INT dim);
	
	
	INT	getN() const;
	INT	getDim() const;
	T	*getP(INT n = 0);
	
	void	clear();
	
	void	setSlave();
	
//------------------------------------------------------------------------

	friend ostream & operator<< <T> (ostream & s, const CIVectors<T> &);	

//------------------------------------------------------------------------
	
private:
INT	dim;					//	dimension of ci-vectors
INT	n;						//	number ci-vectors
SharedMemory	*shared;	//	pointer to shared memory object
T	*p;						//	pointer to data
INT	delFlag;				//	if true deletion neccessary
};


template <class T>	ostream & operator<<(ostream & s, const CIVectors<T> &);	

template <class T>
inline
INT	CIVectors<T>::getN() const
{	return n;	}

template <class T>
inline
INT	CIVectors<T>::getDim() const
{	return dim;	}

template <class T>
inline
T *	CIVectors<T>::getP(INT n)
{	return p+n*dim;	}

template <class T>
inline
T &	CIVectors<T>::operator() (INT i, INT j)
{	return	p[j*dim + i];	}

template <class T>
inline
T &	CIVectors<T>::operator() (INT i)
{	return	p[i];	}



#endif

