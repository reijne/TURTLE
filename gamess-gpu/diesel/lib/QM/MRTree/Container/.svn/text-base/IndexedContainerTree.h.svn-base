//***********************************************************************
//
//	Name:			IndexedContainerTree.h
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


#ifndef __INDEXEDCONTAINERTREE_H
#define __INDEXEDCONTAINERTREE_H

#include "../../../../config.h"


#include "ContainerTree.h"
#include "IndexedContainer.h"

#include "../../Configuration/ConfigurationSAFNr.h"


//template <class ParentType, class ContainedObjectType, class TotalObjectType>
template <class ParentType, class ContainedObjectType>
class IndexedContainerTree : 
	public virtual ContainerTree<ParentType, ContainedObjectType>, 
	public virtual IndexedContainer<ContainedObjectType> {
public:
	IndexedContainerTree() {}
	IndexedContainerTree(INT n) : 
		IndexedContainer<ContainedObjectType>(n) {}	
	IndexedContainerTree(INT n, ParentType * parent) : 
		IndexedContainer<ContainedObjectType>(n),
		Tree<ParentType>(parent)	{}	

private:
};



#endif
