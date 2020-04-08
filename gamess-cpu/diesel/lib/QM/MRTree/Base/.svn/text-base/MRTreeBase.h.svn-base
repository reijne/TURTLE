//***********************************************************************
//
//	Name:			MRTreeBase.h
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


#ifndef __MRTreeBase_h
#define __MRTreeBase_h

#include "../../../../config.h"




#include "../Container/ContainerTree.h"

#include "../../Configuration/ConfigurationSAFNr.h"

struct	MRTreeIterator;

//template <class ParentType, class ContainedObjectType, class TotalObjectType>
template <class ParentType, class ContainedObjectType>
class MRTreeBase : 
	public virtual ContainerTree<ParentType, ContainedObjectType> {
public:
	MRTreeBase() {}
//	MRTreeBase(ParentType parent) : 
//		ContainerTree<ParentType, ContainedObjectType>(parent)	{}
	
	ConfigurationSAFNr<MOType>	getConfigurationSAFNr(INT index) const;



	MRTreeIterator	firstInTree() const;
	void	nextInTree(MRTreeIterator &) const;
	INT	isLastInTree(const MRTreeIterator &) const;
	ConfigurationSAFNr<MOType>	getConfigurationSAFNr(
		const MRTreeIterator & iter) const;

	INT	init(INT, const IrRep *);


private:
};



#endif
