//***********************************************************************
//
//	Name:			BaseContainer.cc
//
//	Description:	
//
//	Author:			Michael Hanrath
//
//	Version:		0.0
//
//	Date:			03.09.1996
//
//
//
//
//
//***********************************************************************


#include <iostream>

#include "BaseContainer.h"

using namespace std;

template <class T>
BaseContainer<T>::BaseContainer(INT _n)
{
	n = _n;
	p = new T*[n];
	memset(p, 0, n*sizeof(T *));
}


template <class T>
BaseContainer<T>::~BaseContainer()
{
	for ( INT i=0 ; i<n ; i++ )
		if ( p[i] )
			delete p[i];

	delete[] p;
}


template <class T>
INT	BaseContainer<T>::getNumberOfLeaves() const
{
INT	leaves = 0;

	for ( INT i=0 ; i<n ; i++ )
		if ( p[i] )
			leaves += p[i]->getNumberOfLeaves();
	return leaves;
}

template <class T>
void	BaseContainer<T>::check()
{
	for ( INT i=0 ; i<n ; i++ )
		if ( p[i] )
		{	cout << "-----------------------" << i << endl;
			p[i]->check();
		}
}

template <class T>
INT	BaseContainer<T>::getNumberOfIntegrals() const
{
INT	integrals = 0;

	for ( INT i=0 ; i<n ; i++ )
		if ( p[i] )
			integrals += p[i]->getNumberOfIntegrals();
	return integrals;
}

template <class T>
INT	BaseContainer<T>::getDataSize() const
{
INT	size = 0;

	for ( INT i=0 ; i<n ; i++ )
		if ( p[i] )
			size += p[i]->getDataSize();
	return size;
}

template <class T>
INT	BaseContainer<T>::allocate()
{
	for ( INT i=0 ; i<n ; i++ )
		if ( p[i] )
			if ( !p[i]->allocate() )
				return 0;
	return 1;
}

template <class T>
void	BaseContainer<T>::setSharedMem(SharedMemory *sharedMem)
{
	for ( INT i=0 ; i<n ; i++ )
		if ( p[i] )
			p[i]->setSharedMem(sharedMem);
}




#include "SymmetryContainer.h"
//#include "NExtContainer.h"
#include "ExtPosContainer.h"
#include "NTupelContainer.h"


template class	BaseContainer<SymmetryContainer>;
//template class	BaseContainer<NExtContainer>;
template class	BaseContainer<ExtPosContainer>;
template class	BaseContainer<NTupelContainer>;
