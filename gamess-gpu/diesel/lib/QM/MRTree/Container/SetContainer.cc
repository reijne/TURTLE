//***********************************************************************
//
//	Name:			SetContainer.cc
//
//	Description:	implements the functionality of a set based
//					on the GNU-C++-Library
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




#include "SetContainer.h"


template <class ContainedObjectType>
SetContainer<ContainedObjectType>::SetContainer(istream &s)
{
//	cout << "SetContainer(istream &s)" << endl;
INT	n;
	s >> n;
//	cout << "n=" << n << endl;
	
	for ( INT i=0 ; i<n ; i++ )
	{
	ContainedObjectType	*p = new ContainedObjectType(s);
//		cout << ";;;;;;;;" << *p << endl;
		add(p);
	}
}


#include "../Sel/TupelStructureSel.h"

#include "../Set/InternalConfsSet.h"
#include "../Set/TupelStructureSet.h"
#include "../Set/extMOsSet.h"
#include "../Base/extEntry.h"



template class SetContainer<TupelStructureSel>;

template class SetContainer<InternalConfsSet>;
template class SetContainer<TupelStructureSet>;
template class SetContainer<extMOsSet>;
template class SetContainer<extEntry>;

