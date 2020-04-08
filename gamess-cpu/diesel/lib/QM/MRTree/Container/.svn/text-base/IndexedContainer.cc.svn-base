//***********************************************************************
//
//	Name:			IndexedContainer.cc
//
//	Description:	
//
//	Author:			Michael Hanrath
//
//	Version:		0.0
//
//	Date:			30. Aug 1998
//
//***********************************************************************




#include "IndexedContainer.h"


#include "../Sel/NExternalsSel.h"
#include "../Diag/NExternalsDiag.h"
#include "../Sel/InternalConfsSel.h"
#include "../Diag/InternalConfsDiag.h"
#include "../Sel/TupelStructureSel.h"
#include "../Diag/TupelStructureDiag.h"
#include "../Sel/extMOsSel.h"
#include "../Diag/extMOsDiag.h"
#include "../Base/extEntry.h"


template <class ContainedObjectType>
IndexedContainer<ContainedObjectType>::IndexedContainer(INT _n)
{
	n = _n;
	p = new ContainedObjectType * [n];	
	memset(p, 0, n*sizeof(ContainedObjectType *));
}

template <class ContainedObjectType>
IndexedContainer<ContainedObjectType>::IndexedContainer(istream &s)
{
//	cout << "IndexedContainer(istream &s)" << endl;
	s >> n;
//	cout << "n=" << n << endl;
	if ( n )
		p = new ContainedObjectType * [n];
	else
		p = NULL;
	
	for ( INT i=0 ; i<n ; i++ )
	{
//		cout << "i=" << i << endl;
		p[i] = new ContainedObjectType(s);
//		cout << "Bi=" << i << endl;
	}
}

template <class ContainedObjectType>
IndexedContainer<ContainedObjectType>::~IndexedContainer()
{
	if ( p )
		delete[] p;
}


//template class IndexedContainer<double>;

template class IndexedContainer<extMOsDiag>;
template class IndexedContainer<TupelStructureDiag>;
template class IndexedContainer<InternalConfsDiag>;


template class IndexedContainer<extMOsSel>;
template class IndexedContainer<TupelStructureSel>;
template class IndexedContainer<InternalConfsSel>;

/*
template <class ContainedObjectType> class ContainerTree;
class test;
template class IndexedContainer<ContainerTree<test *> *>;
*/

