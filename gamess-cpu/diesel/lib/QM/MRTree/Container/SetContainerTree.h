//***********************************************************************
//
//	Name:			SetContainerTree.h
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


#ifndef __SetContainerTree_H
#define __SetContainerTree_H

#include "../../../../config.h"


#include "ContainerTree.h"
#include "SetContainer.h"

#include "../../Configuration/ConfigurationSAFNr.h"


//template <class ParentType, class ContainedObjectType, class TotalObjectType>
template <class ParentType, class ContainedObjectType>
class SetContainerTree : 
	public virtual ContainerTree<ParentType, ContainedObjectType>, 
	public virtual SetContainer<ContainedObjectType> {
public:
	SetContainerTree() {}
	SetContainerTree(ParentType * parent) : 
		Tree<ParentType>(parent)	{}
	virtual ~SetContainerTree();

	// "virtual" constructor:
	virtual SetContainerTree<ParentType, ContainedObjectType> *
		new_SetContainerTree();

	void	copy(SetContainerTree<ParentType, ContainedObjectType> *, INT deep = 0);
	virtual SetContainerTree<ParentType, ContainedObjectType> *
		clone(INT deep = 0);
	

	//----------------------------------------------------------------
	
	SetContainerTree<ParentType, ContainedObjectType> &
		operator |= (SetContainerTree<ParentType, ContainedObjectType>& b);
	SetContainerTree<ParentType, ContainedObjectType> &
		operator &= (SetContainerTree<ParentType, ContainedObjectType>& b);
	SetContainerTree<ParentType, ContainedObjectType> &
		operator -= (SetContainerTree<ParentType, ContainedObjectType>& b);

	//----------------------------------------------------------------

private:
};



template <class ParentType, class ContainedObjectType>
SetContainerTree<ParentType, ContainedObjectType> *
	SetContainerTree<ParentType, ContainedObjectType>::new_SetContainerTree()
{	return new SetContainerTree<ParentType, ContainedObjectType>();	}

template <class ParentType, class ContainedObjectType>
void	SetContainerTree<ParentType, ContainedObjectType>::
	copy(SetContainerTree<ParentType, ContainedObjectType> *c, INT deep)
{
	if ( deep == 0 )
		;
	else
		for ( ContainerIterator i=c->first() ; !c->isLast(i) ; c->next(i) )
			add((*c)[i]);
}

template <class ParentType, class ContainedObjectType>
SetContainerTree<ParentType, ContainedObjectType> *
	SetContainerTree<ParentType, ContainedObjectType>::clone(INT deep)
{	
	SetContainerTree<ParentType, ContainedObjectType> * c = 
		new SetContainerTree<ParentType, ContainedObjectType>();
	c->copy(this, deep);
	return c;
}




#endif
