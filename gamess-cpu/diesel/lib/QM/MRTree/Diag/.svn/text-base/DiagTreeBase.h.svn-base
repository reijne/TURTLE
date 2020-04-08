//***********************************************************************
//
//	Name:			DiagTreeBase.h
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


#ifndef __DiagTreeBase_h
#define __DiagTreeBase_h

#include "../../../../config.h"


#include "../Container/IndexedContainerTree.h"
#include "../Base/MRTreeBase.h"


template <class ParentType, class ContainedObjectType>
class DiagTreeBase :
	public MRTreeBase<ParentType, ContainedObjectType>,
	virtual public IndexedContainerTree<ParentType, ContainedObjectType> {
public:
	DiagTreeBase();
	DiagTreeBase(INT n);
	DiagTreeBase(INT n, ParentType *parent);

	~DiagTreeBase();
	
	
};





#endif
