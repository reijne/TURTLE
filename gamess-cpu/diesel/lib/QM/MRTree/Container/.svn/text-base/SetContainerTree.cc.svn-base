//***********************************************************************
//
//	Name:			SetContainerTree.cc
//
//	Description:	
//
//	Author:			Michael Hanrath
//
//	Version:		0.0
//
//	Date:			21.06.1996
//
//
//
//
//
//***********************************************************************




#include "SetContainerTree.h"


#include "../Sel/NExternalsSel.h"
#include "../Sel/TupelStructureSel.h"

#include "../Set/NExternalsSet.h"
#include "../Set/InternalConfsSet.h"
#include "../Set/TupelStructureSet.h"
#include "../Set/extMOsSet.h"
#include "../Base/extEntry.h"


INT	level = 0;


//	Bug in GNU C++ 2.7.2.3 on AIX: Classes with virtual functions must have
//	at least one virtual function defined non-inlined
template <class ParentType, class ContainedObjectType>
SetContainerTree<ParentType, ContainedObjectType>::~SetContainerTree()
{
}


template <class ParentType, class ContainedObjectType>
SetContainerTree<ParentType, ContainedObjectType> &	
	SetContainerTree<ParentType, ContainedObjectType>::
	operator |= (SetContainerTree<ParentType, ContainedObjectType> & b)
{

//	cout << "start level " << ++level << endl;
//	cout << "A: length=" << length() << endl;
	SetContainer<ContainedObjectType>::operator |= (b);
//	cout << "B: length=" << length() << endl;
	for ( ContainerIterator iter = b.first() ; !b.isLast(iter) ; b.next(iter) )
		if ( b[iter] )
		{
//			cout << "AAAAAA:" << *(operator [] (seek(b[iter]))) << endl;
//			cout << "BBBBBB:" << (*b[iter]) << endl;
			operator [] (seek(b[iter]))->operator |= (*b[iter]);
		}
//	cout << "end level " << level-- << endl;
	return *this;
}


template <class ParentType, class ContainedObjectType>
SetContainerTree<ParentType, ContainedObjectType> &	
	SetContainerTree<ParentType, ContainedObjectType>::
	operator &= (SetContainerTree<ParentType, ContainedObjectType> & b)
{
	SetContainer<ContainedObjectType>::operator &= (b);
	for ( ContainerIterator iter = this->first() ; !this->isLast(iter) ; this->next(iter) )
		if ( b[iter] )
			this->operator [] (iter)->operator &= (*b[seek(this->operator [] (iter))]);
	return *this;
}

#define MEQ \
template <>\
SetContainerTree<ParentType, ContainedObjectType> &	\
	SetContainerTree<ParentType, ContainedObjectType>::\
	operator -= (SetContainerTree<ParentType, ContainedObjectType> & b)\
{\
	SetContainer<ContainedObjectType>::operator -= (b);\
	for ( ContainerIterator iter = b.first() ; !b.isLast(iter) ; b.next(iter) )\
		if ( b[iter] )\
			operator [] (seek(b[iter]))->SetContainer<ContainedContainedObjectType>::operator -= (*b[iter]);\
	return *this;\
}

template<>
SetContainerTree<TupelStructureSet, extEntry> &
	SetContainerTree<TupelStructureSet, extEntry>::
	operator -= (SetContainerTree<TupelStructureSet, extEntry> & b)
{
	SetContainer<extEntry>::operator -= (b);
	for ( ContainerIterator iter = b.first() ; !b.isLast(iter) ; b.next(iter) )
		if ( b[iter] )
			operator [] (seek(b[iter]))->operator -= (*b[iter]);
	return *this;
}



#define ParentType void
#define ContainedObjectType InternalConfsSet
#define ContainedContainedObjectType TupelStructureSet
MEQ
#undef ParentType
#undef ContainedObjectType
#undef ContainedContainedObjectType

#define ParentType InternalConfsSet
#define ContainedObjectType TupelStructureSet
#define ContainedContainedObjectType extMOsSet
MEQ
#undef ParentType
#undef ContainedObjectType
#undef ContainedContainedObjectType

#define ParentType TupelStructureSet
#define ContainedObjectType extMOsSet
#define ContainedContainedObjectType extEntry
MEQ
#undef ParentType
#undef ContainedObjectType
#undef ContainedContainedObjectType


//template class SetContainerTree<NExternalsSel, TupelStructureSel>;

template class SetContainerTree<void, InternalConfsSet>;
template class SetContainerTree<NExternalsSet, TupelStructureSet>;
template class SetContainerTree<InternalConfsSet, extMOsSet>;
template class SetContainerTree<TupelStructureSet, extEntry>;

template class SetContainerTree<TupelStructureBase, extEntry>;
template class SetContainerTree<NExternalsSel, TupelStructureSel>;

