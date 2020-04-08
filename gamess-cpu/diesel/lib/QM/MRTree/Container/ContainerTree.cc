//***********************************************************************
//
//	Name:			ContainerTree.cc
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




#include "ContainerTree.h"


template <class ParentType, class ContainedObjectType>
INT	ContainerTree<ParentType, ContainedObjectType>::getNumberOfLeaves() const
{
INT	numberOfLeaves = 0;
	for ( ContainerIterator iter = this->first() ; !this->isLast(iter) ; this->next(iter) )
	{
//		printf("getNumberOfLeaves() = %d\n",
//			operator [] (iter)->getNumberOfLeaves());
		if ( ((ContainerTree<ParentType, ContainedObjectType> *) this)
				->operator [] (iter) )
			numberOfLeaves += 
				((ContainerTree<ParentType, ContainedObjectType> *) this)
				->operator [] (iter)->getNumberOfLeaves();	
	}
	return numberOfLeaves;
}


template <class ParentType, class ContainedObjectType>
void	ContainerTree<ParentType, ContainedObjectType>::cutTree()
{
	for ( ContainerIterator iter = this->first() ; !this->isLast(iter) ; this->next(iter) )
	{
//		printf("getNumberOfLeaves() = %d\n",
//			operator [] (iter)->getNumberOfLeaves());
		if ( ((ContainerTree<ParentType, ContainedObjectType> *) this)
				->operator [] (iter) )
		{
			if ( ((ContainerTree<ParentType, ContainedObjectType> *) this)
				->operator [] (iter)->getNumberOfLeaves()==0 )
			{
				delete 	((ContainerTree<ParentType, ContainedObjectType> *) this)
				->operator [] (iter);
				((ContainerTree<ParentType, ContainedObjectType> *) this)
				->operator [] (iter) = NULL;
			}
			else
				((ContainerTree<ParentType, ContainedObjectType> *) this)
				->operator [] (iter)->cutTree();	
		}
	}
}


#include "../Diag/TupelStructureDiag.h"
#include "../Diag/InternalConfsDiag.h"
#include "../Diag/NExternalsDiag.h"

#include "../Sel/extMOsSel.h"
#include "../Sel/TupelStructureSel.h"
#include "../Sel/InternalConfsSel.h"
#include "../Sel/NExternalsSel.h"

#include "../Set/extMOsSet.h"
#include "../Set/TupelStructureSet.h"
#include "../Set/InternalConfsSet.h"
#include "../Set/NExternalsSet.h"


template class ContainerTree<InternalConfsDiag, extMOsDiag>;
template class ContainerTree<NExternalsDiag, TupelStructureDiag>;
template class ContainerTree<void, InternalConfsDiag>;

template class ContainerTree<InternalConfsSel, extMOsSel>;
template class ContainerTree<NExternalsSel, TupelStructureSel>;
template class ContainerTree<void, InternalConfsSel>;

template class ContainerTree<InternalConfsSet, extMOsSet>;
template class ContainerTree<NExternalsSet, TupelStructureSet>;
template class ContainerTree<void, InternalConfsSet>;

