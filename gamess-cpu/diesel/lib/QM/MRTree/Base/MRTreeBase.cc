//***********************************************************************
//
//	Name:			MRTreeBase.cc
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




#include "MRTreeBase.h"


#include "../../MRTree/MRTreeIterator.h"

#include "../Sel/extMOsSel.h"
#include "../Sel/TupelStructureSel.h"
#include "../Sel/InternalConfsSel.h"
#include "../Sel/NExternalsSel.h"

#include "../Set/extMOsSet.h"
#include "../Set/TupelStructureSet.h"
#include "../Set/InternalConfsSet.h"
#include "../Set/NExternalsSet.h"

#include "../Diag/extMOsDiag.h"
#include "../Diag/TupelStructureDiag.h"
#include "../Diag/InternalConfsDiag.h"
#include "../Diag/NExternalsDiag.h"



#define InstantiateClassTemplates \
template <>\
MRTreeIterator	MRTreeBase<void, InternalConfs>::\
	firstInTree() const\
{\
MRTreeBase<void, InternalConfs>	*thisNonConst =\
	((MRTreeBase<void, InternalConfs>	*) this);\
	\
MRTreeIterator	mriter;\
	mriter.Citer[0] = thisNonConst->first();\
	\
InternalConfs	*intern = thisNonConst->operator [] (mriter.Citer[0]);\
	mriter.Citer[1] = intern->first();\
\
TupelStructure	*tupel = intern->operator [] (mriter.Citer[1]);\
	mriter.Citer[2] = tupel->first();\
\
	mriter.extMOs = tupel->operator [] (mriter.Citer[2]);\
	mriter.Citer[3] = mriter.extMOs->first();\
	\
	mriter.n = getNumberOfLeaves() - 1;\
	mriter.i = 0;\
	return mriter;\
}\
\
\
template <>\
void	MRTreeBase<void, InternalConfs>::\
	nextInTree(MRTreeIterator &mriter) const\
{\
MRTreeBase<void, InternalConfs>	*thisNonConst =\
	((MRTreeBase<void, InternalConfs>	*) this);\
\
	mriter.n--;\
	mriter.i++;\
	if ( mriter.n<0 )\
		return;\
	mriter.extMOs->next(mriter.Citer[3]);\
	if ( !mriter.extMOs->isLast(mriter.Citer[3]) )\
		return;\
\
	while ( 1 )\
	{\
	InternalConfs	*intern = thisNonConst->operator [] (mriter.Citer[0]);\
	TupelStructure	*tupel = intern->operator [] (mriter.Citer[1]);\
		\
		tupel->next(mriter.Citer[2]);\
		\
		if ( !tupel->isLast(mriter.Citer[2]) )\
		{\
			mriter.extMOs = tupel->operator [] (mriter.Citer[2]);\
			if ( !mriter.extMOs )\
				continue;\
			if ( !mriter.extMOs->getNumberOfElements() )\
				continue;\
			mriter.Citer[3] = mriter.extMOs->first();\
			return;\
		}\
		intern->next(mriter.Citer[1]);\
		if ( !intern->isLast(mriter.Citer[1]) )\
		{\
			tupel = intern->operator [] (mriter.Citer[1]);\
			if ( !tupel )\
				continue;\
			if ( !tupel->getNumberOfElements() )\
				continue;\
			mriter.Citer[2] = tupel->first();\
			mriter.extMOs = tupel->operator [] (mriter.Citer[2]);\
			if ( !mriter.extMOs )\
				continue;\
			if ( !mriter.extMOs->getNumberOfElements() )\
				continue;\
			mriter.Citer[3] = mriter.extMOs->first();\
			return;\
		}\
		do\
		{\
			thisNonConst->next(mriter.Citer[0]);\
			if ( thisNonConst->isLast(mriter.Citer[0]) )\
				break;\
			if ( !thisNonConst->operator [] (mriter.Citer[0]) ) \
				break;\
		} while ( thisNonConst->operator [] (mriter.Citer[0])->getNumberOfElements()==0 );\
		if ( !thisNonConst->isLast(mriter.Citer[0]) )\
		{\
			intern = thisNonConst->operator [] (mriter.Citer[0]);\
			if ( !intern )\
				continue;\
			mriter.Citer[1] = intern->first();\
			tupel = intern->operator [] (mriter.Citer[1]);\
			if ( !tupel )\
				continue;\
			if ( !tupel->getNumberOfElements() )\
				continue;\
			mriter.Citer[2] = tupel->first();\
			mriter.extMOs = tupel->operator [] (mriter.Citer[2]);\
			if ( !mriter.extMOs )\
				continue;\
			if ( !mriter.extMOs->getNumberOfElements() )\
				continue;\
			mriter.Citer[3] = mriter.extMOs->first();\
			return;\
		}\
		cout << "Error in \"void	MRTreeBase<void, InternalConfs>::" << endl;\
		cout << "nextInTree(MRTreeIterator &iter) const\": iteration after tail. " << endl;\
		exit(1);\
	}\
}\


template <class ParentType, class ContainedObjectType>
INT	MRTreeBase<ParentType, ContainedObjectType>::
	isLastInTree(const MRTreeIterator & mriter) const
{
	return mriter.n<0;
}



template <class ParentType, class ContainedObjectType>
ConfigurationSAFNr<MOType>	MRTreeBase<ParentType, ContainedObjectType>::
	getConfigurationSAFNr(const MRTreeIterator & mriter) const
{
/*	return 	operator [] (mriter.Citer[0]) ->
			operator [] (mriter.Citer[1]) ->
			operator [] (mriter.Citer[2]) ->
			getConfigurationSAFNr(mriter.Citer[3]);
*/
	
	return	mriter.extMOs->getConfigurationSAFNr(mriter.Citer[3]);
}



template <class ParentType, class ContainedObjectType>
ConfigurationSAFNr<MOType>	MRTreeBase<ParentType, ContainedObjectType>
	::getConfigurationSAFNr(INT n) const
{
INT	count = 0;
INT	j;
ContainerIterator iter;
MRTreeBase<ParentType, ContainedObjectType>	*thisNonConst =
	((MRTreeBase<ParentType, ContainedObjectType>	*) this);
//	cout << "-------" << endl;
	for ( iter = this->first() ; !this->isLast(iter) ; this->next(iter) )
	{
		j =  thisNonConst->operator [] (iter)->getNumberOfLeaves();
//		cout << "count= " << count << ", j= " << j << "   " << i << endl;
		if ( count+j>n )
			break;
		count += j;
	}
//	return  ((MRTreeBase<ParentType, ContainedObjectType> *)
//			operator [] (iter))->operator() (n-count);
	return  thisNonConst->operator [] (iter)->getConfigurationSAFNr(n-count);
}

template <class ParentType, class ContainedObjectType>
INT	MRTreeBase<ParentType, ContainedObjectType>
	::init(INT n, const IrRep *irrep)
{
	for ( ContainerIterator iter = this->first() ; !this->isLast(iter) ; this->next(iter) )
		if ( this->operator [] (iter) )
			n = this->operator [] (iter)->init(n, irrep);
/*	for ( INT i=0 ; i<getNumberOfElements() ; i++ )
		if ( (*this)[i] )
			n = (*this)[i]->init(n, irrep);
*/	return n;
}


#include "../../Configuration/Configuration.h"

#define	InternalConfs InternalConfsDiag
#define TupelStructure TupelStructureDiag
InstantiateClassTemplates
#undef InternalConfs
#undef TupelStructure

#define	InternalConfs InternalConfsSel
#define TupelStructure TupelStructureSel
InstantiateClassTemplates
#undef InternalConfs
#undef TupelStructure

#define	InternalConfs InternalConfsSet
#define TupelStructure TupelStructureSet
InstantiateClassTemplates
#undef InternalConfs
#undef TupelStructure


template class MRTreeBase
	<InternalConfsSel, extMOsSel>;
template class MRTreeBase
	<NExternalsSel, TupelStructureSel>;
template class MRTreeBase
	<void, InternalConfsSel>;


template class MRTreeBase
	<InternalConfsSet, extMOsSet>;
template class MRTreeBase
	<NExternalsSet, TupelStructureSet>;
template class MRTreeBase
	<void, InternalConfsSet>;


template class MRTreeBase
	<InternalConfsDiag, extMOsDiag>;
template class MRTreeBase
	<NExternalsDiag, TupelStructureDiag>;
template class MRTreeBase
	<void, InternalConfsDiag>;

