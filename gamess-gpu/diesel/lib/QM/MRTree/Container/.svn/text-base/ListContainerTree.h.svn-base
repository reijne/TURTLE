//***********************************************************************
//
//	Name:			ListContainerTree.h
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


#ifndef __ListContainerTree_H
#define __ListContainerTree_H

#include "../../../../config.h"


#include "ContainerTree.h"
#include "ListContainer.h"

#include "../../Configuration/ConfigurationSAFNr.h"


//template <class ParentType, class ContainedObjectType, class TotalObjectType>
template <class ParentType, class ContainedObjectType>
class ListContainerTree : 
	public virtual ContainerTree<ParentType, ContainedObjectType>, 
	public virtual ListContainer<ContainedObjectType> {
public:
	ListContainerTree() {}
	ListContainerTree(ParentType * parent) : 
		Tree<ParentType>(parent)	{}	

private:
};



#endif
