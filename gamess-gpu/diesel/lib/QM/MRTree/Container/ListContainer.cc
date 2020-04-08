//***********************************************************************
//
//	Name:			ListContainer.cc
//
//	Description:	implements the functionality of a singly linked list
//					based on the GNU-C++-Library
//
//	Author:			Michael Hanrath
//
//	Version:		0.0
//
//	Date:			08.03.1997
//
//
//
//
//
//***********************************************************************




#include "ListContainer.h"



/*
template <class ContainedObjectType>
ListContainer<ContainedObjectType>::ListContainer(istream &s) :
  SLList<ContainedObjectType *>
{
//	cout << "ListContainer(istream &s)" << endl;
INT	n;
	s >> n;
//	cout << "n=" << n << endl;
	
	for ( INT i=0 ; i<n ; i++ )
	{
	ContainedObjectType	*p = new ContainedObjectType(s);
		append(p);
	}
}
*/


#include "../Base/extEntry.h"

template class ListContainer<extEntry>;


