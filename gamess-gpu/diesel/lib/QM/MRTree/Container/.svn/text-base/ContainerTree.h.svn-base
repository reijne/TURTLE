//***********************************************************************
//
//	Name:			ContainerTree.h
//
//	Description:	a container tree base class
//
//	Author:			Michael Hanrath
//
//	Version:		0.0
//
//	Date:			14.03.1997
//
//
//
//
//
//***********************************************************************


#ifndef __ContainerTree_H
#define __ContainerTree_H

#include "../../../../config.h"

#include "Tree.h"
#include "Container.h"



template <class ParentType, class ContainedObjectType>
class ContainerTree :
	public virtual Tree<ParentType>,
	public virtual Container<ContainedObjectType> {
public:
	ContainerTree() {}
//	ContainerTree(ParentType parent) : Tree<ParentType>(parent) {}

/*
	// "virtual" constructor:
	virtual ContainerTree<ParentType, ContainedObjectType> *
		new_ContainerTree();

	void	copy(ContainerTree<ParentType, ContainedObjectType> *, INT deep = 0);
	virtual ContainerTree<ParentType, ContainedObjectType> *
		clone(INT deep = 0);
*/	

	INT	getNumberOfLeaves() const;
	void	cutTree();

private:
};




/*
template <class ParentType, class ContainedObjectType>
ContainerTree<ParentType, ContainedObjectType> *
	ContainerTree<ParentType, ContainedObjectType>::new_ContainerTree()
{	return new ContainerTree<ParentType, ContainedObjectType>();	}

template <class ParentType, class ContainedObjectType>
void	ContainerTree<ParentType, ContainedObjectType>::
	copy(ContainerTree<ParentType, ContainedObjectType> *c, INT deep)
{
	if ( deep == 0 )
		;
	else
	{
		for ( ContainerIterator i=c->first() ; !c->isLast(i) ; c->next(i) )
//			append((*c)[i]);
			;
	}
}

template <class ParentType, class ContainedObjectType>
ContainerTree<ParentType, ContainedObjectType> *
	ContainerTree<ParentType, ContainedObjectType>::clone(INT deep)
{	
	ContainerTree<ParentType, ContainedObjectType> * c = 
		new ContainerTree<ParentType, ContainedObjectType>();
	c->copy(this, deep);
	return c;
}
*/
#endif
