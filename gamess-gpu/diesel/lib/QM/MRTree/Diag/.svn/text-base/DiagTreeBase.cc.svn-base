//***********************************************************************
//
//	Name:			DiagTreeBase.cc
//
//	Description:	
//
//	Author:			Michael Hanrath
//
//	Version:		0.0
//
//	Date:			15.03.1997
//
//
//
//
//
//***********************************************************************




#include "DiagTreeBase.h"





template <class ParentType, class ContainedObjectType>
DiagTreeBase<ParentType, ContainedObjectType>::DiagTreeBase() 
{
}

template <class ParentType, class ContainedObjectType>
DiagTreeBase<ParentType, ContainedObjectType>::DiagTreeBase(INT n) : 
		IndexedContainerTree<ParentType, ContainedObjectType>(n)
{
}

template <class ParentType, class ContainedObjectType>
DiagTreeBase<ParentType, ContainedObjectType>::DiagTreeBase(
	INT n, ParentType *parent) :
		Tree<ParentType>(parent),
		IndexedContainer<ContainedObjectType>(n)
{
}

template <class ParentType, class ContainedObjectType>
DiagTreeBase<ParentType, ContainedObjectType>::~DiagTreeBase()
{
	for ( ContainerIterator iter = this->first() ; !this->isLast(iter) ; this->next(iter) )
		if ( this->operator [] (iter) )
			delete this->operator [] (iter);
}


#include "../../Configuration/Configuration.h"
#include "extMOsDiag.h"
#include "TupelStructureDiag.h"
#include "InternalConfsDiag.h"


class NExternalsDiag;
template class DiagTreeBase
	<InternalConfsDiag, extMOsDiag>;
template class DiagTreeBase
	<NExternalsDiag, TupelStructureDiag>;
template class DiagTreeBase
	<void, InternalConfsDiag>;

